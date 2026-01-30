function [h_wo,h_w] = plot_freqresp_eqchange(D_wo, D_w, win, param, figTitle, coeffs_ch)
% plot_freqresp_eqchange
%   等化前(woTHP)と等化後(wTHP)で「有効チャネル」をLS推定し、
%   周波数応答（Magnitude + Phase）を比較表示する。
%
% 入力:
%   D_wo, D_w : evaluate_from_io の第2出力（D.tx_al, D.rx_soft_al を含む）
%   win       : first_C, C_max（推定したいtap範囲）
%   param     : Tc_but, cutoff_coeff
%   figTitle  : 図のタイトル文字列
%   coeffs_ch : (任意) 物理チャネル係数を重ねるなら渡す（coeffs.ChannelCoeff）

mRange = win.first_C:win.C_max;

% --- 有効チャネル(woTHP / wTHP) をLS推定 ---
h_wo = estimate_fir_ls_sym(D_wo.tx_al, D_wo.rx_soft_al, mRange);
h_w  = estimate_fir_ls_sym(D_w.tx_al,  D_w.rx_soft_al,  mRange);

% m=0 で正規化して形を比較（主タップ=1）
idx0 = 1 - mRange(1);  % m=0 が入っている位置
h_wo = h_wo / h_wo(idx0);
h_w  = h_w  / h_w(idx0);

% --- 周波数軸（シンボルレート基準） ---
Tc = param.Tc_but / param.cutoff_coeff;
fs_sym = 1/Tc;

Nfft = 4096;
f = (-Nfft/2:Nfft/2-1)/Nfft * fs_sym;

H_wo = fftshift(fft(h_wo, Nfft));
H_w  = fftshift(fft(h_w,  Nfft));

% (任意) 物理チャネル
has_ch = (nargin >= 6 && ~isempty(coeffs_ch));

if has_ch
    g = coeffs_ch.ChannelCoeff(:).';
    if numel(g) ~= numel(mRange)
        has_ch = false;
    else
        g  = g / g(idx0);
        Hg = fftshift(fft(g, Nfft));
    end
end

% =========================================================
% (1) Magnitude [dB]
% =========================================================
figure; hold on;
plot(f/1e9, 20*log10(abs(H_wo)+1e-12), 'DisplayName','effective woTHP (LS)');
plot(f/1e9, 20*log10(abs(H_w )+1e-12), 'DisplayName','effective wTHP  (LS)');
if has_ch
    plot(f/1e9, 20*log10(abs(Hg)+1e-12), 'DisplayName','ChannelCoeff (given)');
end
grid on;
xlabel('Frequency [GHz]'); ylabel('Magnitude [dB]');
title(['Freq response MAG (before/after eq): ' char(figTitle)]);
legend('show','Location','best');

% =========================================================
% (2) Phase [rad]  (unwrap)
% =========================================================
figure; hold on;
plot(f/1e9, unwrap(angle(H_wo)), 'DisplayName','phase woTHP (LS)');
plot(f/1e9, unwrap(angle(H_w )), 'DisplayName','phase wTHP  (LS)');
if has_ch
    plot(f/1e9, unwrap(angle(Hg)), 'DisplayName','phase ChannelCoeff (given)');
end
grid on;
xlabel('Frequency [GHz]'); ylabel('Phase [rad] (unwrap)');
title(['Freq response PHASE (before/after eq): ' char(figTitle)]);
legend('show','Location','best');

% =========================================================
% (3) 参考：Group delay [s]（位相の傾き）
%     角周波数 w=2*pi*f なので tau = - d(phi)/dw
% =========================================================
w = 2*pi*f;  % [rad/s]
phi_wo = unwrap(angle(H_wo));
phi_w  = unwrap(angle(H_w));

tau_wo = -gradient(phi_wo, w);  % [s]
tau_w  = -gradient(phi_w , w);  % [s]

figure; hold on;
plot(f/1e9, tau_wo*1e12, 'DisplayName','group delay woTHP (LS)');
plot(f/1e9, tau_w *1e12, 'DisplayName','group delay wTHP  (LS)');
if has_ch
    phi_g = unwrap(angle(Hg));
    tau_g = -gradient(phi_g, w);
    plot(f/1e9, tau_g*1e12, 'DisplayName','group delay ChannelCoeff (given)');
end
grid on;
xlabel('Frequency [GHz]'); ylabel('Group delay [ps]');
title(['Group delay (from phase slope): ' char(figTitle)]);
legend('show','Location','best');

% =========================================================
% (4) time-domain taps
% =========================================================
figure; hold on;
stem(mRange , real(h_wo), 'DisplayName','h\_eff woTHP (real)');
stem(mRange + param.first_F, real(h_w ), 'DisplayName','h\_eff wTHP  (real)');
grid on; xlabel('tap index m'); ylabel('tap value');
title(['Effective taps (LS) [Real]: ' char(figTitle)]);
legend('show','Location','best');

if ~isreal(h_wo) || ~isreal(h_w)
    figure; hold on;
    stem(mRange , imag(h_wo), 'DisplayName','h\_eff woTHP (imag)');
    stem(mRange + param.first_F, imag(h_w ), 'DisplayName','h\_eff wTHP  (imag)');
    grid on; xlabel('tap index m'); ylabel('tap value');
    title(['Effective taps (LS) [Imag]: ' char(figTitle)]);
    legend('show','Location','best');
end

print_metrics = @(name, hvec) local_print_metrics(name, hvec, mRange);

print_metrics("h_wo (LS)", h_wo);
print_metrics("h_w  (LS)", h_w);

if has_ch
    print_metrics("ChannelCoeff (given)", g);
end
end

function local_print_metrics(name, h, m)
    H0  = sum(h);
    Hpi = sum(h .* ((-1).^m));
    Se  = sum(h(mod(m,2)==0));
    So  = sum(h(mod(m,2)~=0));
    fprintf('[%s] H(0)=%.4f, |H(pi)|=%.4f, ratio(pi/0)=%.2f dB\n', ...
        name, real(H0), abs(Hpi), 20*log10(abs(Hpi)/abs(H0)));
    fprintf('      Se=%.4f, So=%.4f  => H0=Se+So=%.4f, Hpi=Se-So=%.4f\n', ...
        real(Se), real(So), real(Se+So), real(Se-So));
end
