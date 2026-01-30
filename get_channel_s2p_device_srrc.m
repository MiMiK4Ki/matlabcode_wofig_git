function [ChanFrame, info, param_out] = get_channel_s2p_device_srrc(param_in, win, s2pFile, Zs, Zl, opt)
% get_channel_s2p_device_srrc
%   A1(get_channel_synthetic) と同じ規約で、
%     - s2p -> H_dev(f) -> h_dev(t)（dt刻みの連続近似）
%     - sRRC -> comm.RaisedCosineTransmitFilter で h_srrc(t)（dt刻み）
%     - 連続畳み込み: h_total = conv(h_dev, h_srrc) * dt
%   を作って ChanFrame.h_rx に入れる。
%
% 仕様:
%   - Bdev = -3dB帯域（低域基準から -3dB）
%   - Bsrrc = Bdev * cutoff_coeff
%   - Tc = 1/(2*Bsrrc)
%   - dt = Tc/NoSpS
%   - sRRC は通過域(DC)が 1 になるように integral 正規化

% =========================
% (0) s2p load & H_meas(f)
% =========================
sp = sparameters(s2pFile);                 % RF Toolbox
f_meas = sp.Frequencies(:);                % [Hz]
if any(diff(f_meas) <= 0)
    error('s2p frequency must be strictly increasing.');
end

Htmp = s2tf(sp, Zs, Zl);                   % voltage transfer
H_meas = squeeze(Htmp);                    % ensure vector
H_meas = H_meas(:);

f_start = f_meas(1);
f_end   = f_meas(end);

% =========================
% (1) デバイス -3dB 帯域 Bdev
% =========================
Bdev = estimate_bw_3dB_from_lowfreq(f_meas, H_meas, opt.Kfit);

% Bsrrc = Bdev * cutoff
Bsrrc = Bdev * param_in.cutoff_coeff;

% Tc = 1/(2*Bsrrc) を強制
Tc = 1/(2*Bsrrc);

% dt = Tc/NoSpS
NoSpS = param_in.NoSpS;
dt = Tc / NoSpS;
fs = 1/dt;

% 既存コード互換: Tc = Tc_but/cutoff を満たすよう Tc_but を更新
% Tc_but = Tc * cutoff = 1/(2*Bdev)
param_out = param_in;
param_out.Tc_but = 1/(2*Bdev);
param_out.f_butterworth_cutoff = (1/param_out.Tc_but)/2;

% Nyquist check: s2p測定帯域が fs/2 を超えると IFFTが破綻
if f_end > fs/2
    error(['Nyquist violation: f_end=%.3f GHz > fs/2=%.3f GHz.\n' ...
        'Fix: increase NoSpS (higher fs) or reduce Bsrrc definition.'], ...
        f_end/1e9, (fs/2)/1e9);
end

% =========================
% (2) sRRC impulse (dt刻み連続近似) by comm.RaisedCosineTransmitFilter
% =========================
hFLpS = param_in.hFLpS;            % half span [symbols]
spanSym = 2*hFLpS;                 % full span [symbols]
beta = param_in.beta;

rcTx = comm.RaisedCosineTransmitFilter( ...
    'RolloffFactor', beta, ...
    'FilterSpanInSymbols', spanSym, ...
    'OutputSamplesPerSymbol', NoSpS, ...
    'Gain', 1);



reset(rcTx);

% impulse symbol + zeros to flush
x_imp = [1; zeros(spanSym,1)];     % length spanSym+1 symbols
y_imp = rcTx(x_imp);              % length (spanSym+1)*NoSpS samples

Lsrrc = spanSym*NoSpS + 1;         % desired impulse length
h_srrc = y_imp(1:Lsrrc);           % causal indexing, peak at center
% figure;
% plot(abs(fft(h_srrc)))
% figure;
% plot(y_imp)
% A1相当の時間軸（参考）
t_rrc = (-hFLpS*Tc : dt : hFLpS*Tc); %#ok<NASGU>

% passband=1: DCゲインを 1 に正規化（連続の H(0)=∫h dt）
H0_srrc = sum(h_srrc) * dt;
if abs(H0_srrc) < 1e-15
    error('sRRC DC gain is ~0 (sum(h)*dt). Check parameters.');
end
h_srrc = h_srrc / H0_srrc;         % so that H_srrc(0)=1
% figure;
% plot(abs(fft(h_srrc)).*dt)
% =========================
% (3) H_dev(f) を dt に整合する一様グリッドへ -> IFFTで h_dev(t)
% =========================
df_meas = min(diff(f_meas));
Nfft_min = ceil( 1/(df_meas*dt) ); % df = 1/(Nfft*dt) <= df_meas にしたい
if opt.padPow2
    Nfft = 2^nextpow2(Nfft_min);
else
    Nfft = Nfft_min;
end
if mod(Nfft,2) ~= 0
    Nfft = Nfft + 1;              % 偶数にする
end

df = 1/(Nfft*dt);
Npos = Nfft/2 + 1;
f_pos = (0:Npos-1).' * df;        % [0..fs/2]

idx_fill = (f_pos < f_start);
idx_meas = (f_pos >= f_start) & (f_pos <= f_end);
idx_high = (f_pos > f_end);

% measured interpolation
H_pos_interp = interp1(f_meas, H_meas, f_pos, "linear", "extrap");

% low-frequency fill (optional)
H_pos = zeros(size(f_pos));
if opt.doFillLF
    H_fill = lowfreq_fill_model(f_meas, H_meas, f_pos(idx_fill), opt.Kfit);
    H_pos(idx_fill) = H_fill;
else
    % simplest: hold first measured value
    H_pos(idx_fill) = H_meas(1);
end

H_pos(idx_meas) = H_pos_interp(idx_meas);
H_pos(idx_high) = 0;

% figure
% plot(real(H_pos))
% hold on
% plot(imag(H_pos))

% enforce DC real
H_pos(1) = real(H_pos(1));

% conjugate-symmetric spectrum -> real impulse
H_full = [H_pos; conj(H_pos(end-1:-1:2))];   % length Nfft

% continuous-time scaling: h(t) ≈ Nfft*df * ifft(H)
h_dev_full = (Nfft*df) * ifft(H_full, "symmetric");
h_dev_full = h_dev_full(:);

% =========================
% (4) デバイスIRを A1 と同じ「±hFLpS*Tc 相当」に切り出す（軽量化&整合）
% =========================
[~, idx_peak] = max(abs(h_dev_full));
halfSamps = hFLpS*NoSpS;

i1 = idx_peak - halfSamps;
i2 = idx_peak + halfSamps;

padL = 0; padR = 0;
if i1 < 1
    padL = 1 - i1; i1 = 1;
end
if i2 > numel(h_dev_full)
    padR = i2 - numel(h_dev_full); i2 = numel(h_dev_full);
end

h_dev = [zeros(padL,1); h_dev_full(i1:i2); zeros(padR,1)];
h_dev = h_dev(:).';             % row

% =========================
% (5) 連続畳み込みで合成（A1と同じ conv*dt）
% =========================
h_total = conv(h_dev, h_srrc(:).') * dt;

% =========================
% (6) ChanFrame
% =========================
ChanFrame = struct();
ChanFrame.h_rx  = h_total;      % *dt 済み
ChanFrame.dt    = dt;
ChanFrame.pulse = h_srrc(:).';  % sRRC（*dt無しの連続近似サンプル）
ChanFrame.h_dev = h_dev;        % device impulse (cut)
ChanFrame.meta  = struct('f_meas',f_meas,'H_meas',H_meas,'df',df,'Nfft',Nfft);

info = struct();
info.Bdev_Hz  = Bdev;
info.Bsrrc_Hz = Bsrrc;
info.Tc       = Tc;
info.dt       = dt;
info.fs       = fs;
info.df       = df;
info.Nfft     = Nfft;
info.f_start  = f_start;
info.f_end    = f_end;

% =========================
% (7) plot (optional, no subplot)
% =========================
if opt.doPlot
    % magnitude
    figure; hold on; grid on;
    plot(f_meas/1e9, 20*log10(max(abs(H_meas),1e-15)), '--', 'DisplayName','Measured H(f)');
    plot(f_pos(idx_meas)/1e9, 20*log10(max(abs(H_pos(idx_meas)),1e-15)), '-', 'LineWidth',1.2, ...
        'DisplayName','Used (interp)');
    if any(idx_fill)
        plot(f_pos(idx_fill)/1e9, 20*log10(max(abs(H_pos(idx_fill)),1e-15)), '-', 'LineWidth',1.8, ...
            'DisplayName','Filled (LF)');
    end
    xline(Bdev/1e9,'--','B_{dev}(-3dB)','LabelVerticalAlignment','bottom');
    xlabel('f [GHz]'); ylabel('|H_{dev}(f)| [dB]');
    title('Device transfer from S2P (with LF fill highlighted)');
    legend('Location','best');

    % impulse (device cut)
    % time axes (device, sRRC)
    t_dev  = (-hFLpS*Tc : dt : hFLpS*Tc);
    t_srrc = (-hFLpS*Tc : dt : hFLpS*Tc);

    % ---- total time axis (IMPORTANT) ----
    % If h_total = conv(h_dev, h_srrc)*dt  (continuous convolution approximation)
    t0_dev  = t_dev(1);
    t0_srrc = t_srrc(1);

    t_total = (t0_dev + t0_srrc) + (0:numel(h_total)-1)*dt;

    % ---- plot ----
    figure;
    plot(t_dev*1e9,   real(h_dev),   'LineWidth', 1.2); grid on; hold on;
    plot(t_srrc*1e9,  real(h_srrc),  'LineWidth', 1.2);
    plot(t_total*1e9, real(h_total), 'LineWidth', 1.2);

    xlabel('t [ns]'); ylabel('h(t)');
    legend('device','srrc','total');
    title('device, sRRC, and total impulse responses');

    % =========================
    % Phase response (device) 追加
    % =========================
    figure; hold on; grid on;

    % measured phase
    phi_meas = unwrap(angle(H_meas));

    % constructed H_pos phase (0..fs/2 on uniform grid)
    phi_pos  = unwrap(angle(H_pos));

    % plot: measured (only where measured exists)
    plot(f_meas/1e9, phi_meas, '--', 'DisplayName','Measured phase');

    % plot: used(interp) region
    plot(f_pos(idx_meas)/1e9, phi_pos(idx_meas), '-', 'LineWidth', 1.2, ...
        'DisplayName','Used (interp) phase');

    % plot: filled LF region (if enabled)
    if any(idx_fill)
        plot(f_pos(idx_fill)/1e9, phi_pos(idx_fill), '-', 'LineWidth', 1.8, ...
            'DisplayName','Filled (LF) phase');
    end

    % helpful vertical lines
    xline(f_start/1e9,'--','f_{start}','LabelVerticalAlignment','bottom');
    xline(Bdev/1e9,'--','B_{dev}(-3dB)','LabelVerticalAlignment','bottom');

    xlabel('f [GHz]'); ylabel('Phase [rad]');
    title('Device phase response from S2P (with LF fill highlighted)');
    legend('Location','best');
end
end

function Bdev = estimate_bw_3dB_from_lowfreq(f, H, Kfit)
% 低域基準(ref)から -3dB になる周波数を Bdev とする
K = min(max(5, Kfit), numel(f));
magdB = 20*log10(max(abs(H), 1e-15));
refdB = mean(magdB(1:K));
thr   = refdB - 3;

idx = find(magdB <= thr, 1, 'first');
if isempty(idx)
    Bdev = f(end);
else
    Bdev = f(idx);
end
end

function H_fill = lowfreq_fill_model(f_meas, H_meas, f_fill, Kfit)
% あなたが提示した「低域で遅延推定＋log|H|をf^2で近似」の簡略版
% （補完は 0<=f<f_start にのみ使う）
if isempty(f_fill)
    H_fill = complex([]);
    return;
end

K = min(max(5, Kfit), numel(f_meas));
fK = f_meas(1:K);
HK = H_meas(1:K);

% delay estimate from phase slope
phiK = unwrap(angle(HK));
p = polyfit(fK, phiK, 1);
tau = -p(1)/(2*pi);

% de-embed delay
H_de_K = HK .* exp(1j*2*pi*fK*tau);

% choose DC phase (0 or pi) to make it (mostly) real
phi0_est = angle(mean(H_de_K));
phi0 = 0;
if abs(wrapToPi(phi0_est)) > (pi/2)
    phi0 = pi;
end

% smooth magnitude near DC: log|H_de| vs f^2
x = (fK.^2);
y = log(max(abs(H_de_K), 1e-300));
c = polyfit(x, y, 1);

mag_fill = exp(polyval(c, f_fill.^2));
H_de_fill = mag_fill .* exp(1j*phi0);
H_fill = H_de_fill .* exp(-1j*2*pi*f_fill*tau);

% DC real is enforced outside
end
