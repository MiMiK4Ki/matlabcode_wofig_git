%% ================================================================
% s2p から:
%   (1) スライドの式: S -> ABCD -> (V2/V1, V2/Vs)
%   (2) MATLAB s2tf (option 1/2)
% を重ねて「同一か？」を確認し、
% さらに 0Hz までローパス整合に埋めて実インパルス応答を作る
%
% 重要:
%  - subplot は使わない (figure を分ける)
%  - 周波数応答で「補正(=0〜f_startを生成)した箇所」を色分け
%% ================================================================

clear; close all; clc;

%% ---- User settings ----
s2pFile = "Filter response.s2p";

% Source / load impedances you want for voltage-transfer definitions
Zs = 50;     % [ohm]
Zl = 50;     % [ohm]

% Low-frequency fill model settings
Kfit      = 25;   % use first K points to estimate delay + DC trend
polyOrder = 1;    % fit log|H_de| vs f^2 near DC

% Uniform grid / IFFT settings
padFactor    = 8;         % >=1
interpMethod = "linear";  % "linear" is safest for complex data

%% ---- Load s2p ----
sp = sparameters(s2pFile);
f_meas = sp.Frequencies(:);   % [Hz]
if any(diff(f_meas) <= 0)
    error("Frequency vector must be strictly increasing.");
end

S = sp.Parameters;            % 2x2xN complex
S11 = squeeze(S(1,1,:));
S12 = squeeze(S(1,2,:));
S21 = squeeze(S(2,1,:));
S22 = squeeze(S(2,2,:));

% Reference impedance of the S-parameters (Touchstone Z0)
Z0ref = sp.Impedance;
if isscalar(Z0ref)
    Z0 = Z0ref * ones(size(f_meas));
else
    % If Z0 is frequency-dependent (rare in s2p), make it a column vector
    Z0 = Z0ref(:);
    if numel(Z0) ~= numel(f_meas)
        error("Unsupported sp.Impedance size. Expect scalar or length(f_meas).");
    end
end

f_start = f_meas(1);
fmax    = f_meas(end);

%% ================================================================
% (A) Slide formula: S -> ABCD (manual) -> voltage transfer
%     ABCD formula used (standard, same as your slide):
%     A = ((1+S11)(1-S22)+S12S21)/(2S21)
%     B = Z0((1+S11)(1+S22)-S12S21)/(2S21)
%     C = (1/Z0)((1-S11)(1-S22)-S12S21)/(2S21)
%     D = ((1-S11)(1+S22)+S12S21)/(2S21)
%% ================================================================

A = ((1+S11).*(1-S22) + (S12.*S21)) ./ (2.*S21);
B = (Z0) .* (((1+S11).*(1+S22) - (S12.*S21)) ./ (2.*S21));
C = (1./Z0) .* (((1-S11).*(1-S22) - (S12.*S21)) ./ (2.*S21));
D = ((1-S11).*(1+S22) + (S12.*S21)) ./ (2.*S21);

% Slide (DSP definition): V2/V1 = ZL/(A*ZL + B)
H_V2V1 = Zl ./ (A.*Zl + B);

% Slide (Cadence definition): V2/Vs = ZL/(A*ZL + B + C*Zs*ZL + D*Zs)
H_V2Vs = Zl ./ (A.*Zl + B + C.*Zs.*Zl + D.*Zs);

%% ================================================================
% (B) MATLAB built-ins
%% ================================================================
% s2tf options (MathWorks):
%   option=1 (default): gain from incident voltage Va -> output voltage
%   option=2           : gain from source voltage Vs -> output voltage
H_s2tf_opt1 = s2tf(sp, Zs, Zl, 1);
H_s2tf_opt2 = s2tf(sp, Zs, Zl, 2);

% ABCD check: compare manual ABCD vs MATLAB s2abcd
% s2abcd wants s_params array and z0 (scalar or vector)
try
    ABCD_mat = s2abcd(S, Z0);             % supports vector length(f_meas) in newer versions
catch
    ABCD_mat = s2abcd(S, mean(Z0));       % fallback: scalar Z0
end
A_m = squeeze(ABCD_mat(1,1,:));
B_m = squeeze(ABCD_mat(1,2,:));
C_m = squeeze(ABCD_mat(2,1,:));
D_m = squeeze(ABCD_mat(2,2,:));

% Print ABCD consistency check
fprintf("=== ABCD manual vs s2abcd check ===\n");
fprintf("max|A-A_m| = %.3e\n", max(abs(A - A_m)));
fprintf("max|B-B_m| = %.3e\n", max(abs(B - B_m)));
fprintf("max|C-C_m| = %.3e\n", max(abs(C - C_m)));
fprintf("max|D-D_m| = %.3e\n\n", max(abs(D - D_m)));

%% ================================================================
% (C) Compare transfer functions at measured frequencies (NO fill yet)
%% ================================================================
fprintf("=== Transfer function comparison (measured frequency points) ===\n");
rel = @(x,y) max(abs(x-y) ./ max(abs(y), 1e-12));
fprintf("relErr( H_V2Vs  vs s2tf option2 ) = %.3e  (期待: 小)\n", rel(H_V2Vs, H_s2tf_opt2));
fprintf("relErr( H_V2V1  vs s2tf option2 ) = %.3e  (一般に一致しない)\n", rel(H_V2V1, H_s2tf_opt2));
fprintf("relErr( H_V2V1  vs s2tf option1 ) = %.3e  (一般に一致しない)\n\n", rel(H_V2V1, H_s2tf_opt1));

% ---- Plot: magnitude (measured) ----
figure;
hold on; grid on;
plot(f_meas/1e9, 20*log10(max(abs(H_V2V1),1e-15)), 'LineWidth', 1.2, ...
    'DisplayName','ABCD: V_2/V_1 (slide DSP definition)');
plot(f_meas/1e9, 20*log10(max(abs(H_V2Vs),1e-15)), 'LineWidth', 1.2, ...
    'DisplayName','ABCD: V_2/V_s (slide Cadence definition)');
plot(f_meas/1e9, 20*log10(max(abs(H_s2tf_opt1),1e-15)), '--', 'LineWidth', 1.2, ...
    'DisplayName','s2tf option 1 (V_a \rightarrow V_L)');
plot(f_meas/1e9, 20*log10(max(abs(H_s2tf_opt2),1e-15)), '--', 'LineWidth', 1.2, ...
    'DisplayName','s2tf option 2 (V_s \rightarrow V_L)');
xlabel('f [GHz]');
ylabel('|H(f)| [dB]');
title('Transfer functions at measured frequencies');
legend('Location','best');

% ---- Plot: phase (measured) ----
figure;
hold on; grid on;
plot(f_meas/1e9, unwrap(angle(H_V2V1)), 'LineWidth', 1.2, ...
    'DisplayName','ABCD: V_2/V_1');
plot(f_meas/1e9, unwrap(angle(H_V2Vs)), 'LineWidth', 1.2, ...
    'DisplayName','ABCD: V_2/V_s');
plot(f_meas/1e9, unwrap(angle(H_s2tf_opt1)), '--', 'LineWidth', 1.2, ...
    'DisplayName','s2tf opt1');
plot(f_meas/1e9, unwrap(angle(H_s2tf_opt2)), '--', 'LineWidth', 1.2, ...
    'DisplayName','s2tf opt2');
xlabel('f [GHz]');
ylabel('Phase [rad]');
title('Phase at measured frequencies');
legend('Location','best');

%% ================================================================
% (D) Build real impulse responses:
%     - make uniform grid including DC
%     - fill 0..f_start using lowpass-consistent model
%     - create conjugate-symmetric spectrum and IFFT
%% ================================================================

[f_pos, Hpos_V2V1, idx_fill, idx_meas, t, h_V2V1, df, tau1] = ...
    buildRealImpulseFromMeasuredH(f_meas, H_V2V1, Kfit, polyOrder, padFactor, interpMethod);

[~,     Hpos_s2tf2, ~,      ~,      ~, h_s2tf2, ~,  tau2] = ...
    buildRealImpulseFromMeasuredH(f_meas, H_s2tf_opt2, Kfit, polyOrder, padFactor, interpMethod);

fprintf("=== IFFT grid info ===\n");
fprintf("f_start = %.6f GHz\n", f_start/1e9);
fprintf("df      = %.6f MHz\n", df/1e6);
fprintf("Npos    = %d  (0..%.6f GHz)\n", numel(f_pos), f_pos(end)/1e9);
fprintf("dt      = %.6f ps\n", (t(2)-t(1))*1e12);
fprintf("Estimated delay tau (V2/V1)  = %.6f ns\n", tau1*1e9);
fprintf("Estimated delay tau (s2tf2)  = %.6f ns\n\n", tau2*1e9);

%% ---- Plot magnitude with corrected region highlighted (V2/V1) ----
figure;
hold on; grid on;
plot(f_meas/1e9, 20*log10(max(abs(H_V2V1),1e-15)), '--', 'LineWidth', 1.2, ...
    'DisplayName','Measured samples (ABCD V_2/V_1)');

plot(f_pos(idx_meas)/1e9, 20*log10(max(abs(Hpos_V2V1(idx_meas)),1e-15)), '.', ...
    'DisplayName','Used from measured interpolation (f \ge f_{start})');

plot(f_pos(idx_fill)/1e9, 20*log10(max(abs(Hpos_V2V1(idx_fill)),1e-15)), '-', 'LineWidth', 1.8, ...
    'DisplayName','Filled/extrapolated (0 \le f < f_{start})');

xline(f_start/1e9, '--', 'f_{start}', 'LabelVerticalAlignment','bottom');
xlabel('f [GHz]'); ylabel('|H(f)| [dB]');
title('ABCD V_2/V_1 : magnitude (filled region highlighted)');
legend('Location','best');

%% ---- Plot phase with corrected region highlighted (V2/V1) ----
figure;
hold on; grid on;
phi_pos = unwrap(angle(Hpos_V2V1));
plot(f_pos(idx_meas)/1e9, phi_pos(idx_meas), '.', 'DisplayName','Interp region');
plot(f_pos(idx_fill)/1e9, phi_pos(idx_fill), '-', 'LineWidth', 1.8, 'DisplayName','Filled region');
xline(f_start/1e9, '--', 'f_{start}', 'LabelVerticalAlignment','bottom');
xlabel('f [GHz]'); ylabel('Phase [rad]');
title('ABCD V_2/V_1 : phase (filled region highlighted)');
legend('Location','best');

%% ---- Plot magnitude with corrected region highlighted (s2tf option2) ----
figure;
hold on; grid on;
plot(f_meas/1e9, 20*log10(max(abs(H_s2tf_opt2),1e-15)), '--', 'LineWidth', 1.2, ...
    'DisplayName','Measured samples (s2tf opt2)');

plot(f_pos(idx_meas)/1e9, 20*log10(max(abs(Hpos_s2tf2(idx_meas)),1e-15)), '.', ...
    'DisplayName','Used from measured interpolation (f \ge f_{start})');

plot(f_pos(idx_fill)/1e9, 20*log10(max(abs(Hpos_s2tf2(idx_fill)),1e-15)), '-', 'LineWidth', 1.8, ...
    'DisplayName','Filled/extrapolated (0 \le f < f_{start})');

xline(f_start/1e9, '--', 'f_{start}', 'LabelVerticalAlignment','bottom');
xlabel('f [GHz]'); ylabel('|H(f)| [dB]');
title('s2tf option2 : magnitude (filled region highlighted)');
legend('Location','best');

%% ---- Plot impulse responses (overlay) ----
figure;
hold on; grid on;
plot(t*1e9, h_V2V1, 'LineWidth', 1.2, 'DisplayName','h(t) from ABCD V_2/V_1');
plot(t*1e9, h_s2tf2, 'LineWidth', 1.2, 'DisplayName','h(t) from s2tf option2');
xlabel('t [ns]');
ylabel('h(t)');
title('Real impulse responses (conjugate-symmetric IFFT)');
legend('Location','best');

%% ================================================================
% Local function(s)
%% ================================================================

function [f_pos, H_pos, idx_fill, idx_meas, t, h, df, tau] = ...
    buildRealImpulseFromMeasuredH(f_meas, H_meas, Kfit, polyOrder, padFactor, interpMethod)

    f_meas = f_meas(:);
    H_meas = H_meas(:);

    f_start = f_meas(1);
    fmax    = f_meas(end);

    % ---- choose df (uniform grid) ----
    df = min(diff(f_meas));
    if ~isfinite(df) || df <= 0
        error("Invalid df.");
    end

    Npos = floor(fmax/df) + 1;
    f_pos = (0:(Npos-1)).' * df;

    idx_fill = (f_pos < f_start);
    idx_meas = ~idx_fill;

    % ---- interpolate measured region only (NO extrap) ----
    H_pos_interp = interp1(f_meas, H_meas, f_pos, interpMethod);

    % ---- low-frequency fill (0..f_start) ----
    K = min(Kfit, numel(f_meas));
    fK = f_meas(1:K);
    HK = H_meas(1:K);

    % (1) estimate delay tau from phase slope
    phiK = unwrap(angle(HK));
    p = polyfit(fK, phiK, 1);      % phi ≈ p(1)*f + p(2)
    tau = -p(1)/(2*pi);

    % (2) de-embed delay
    H_de_K = HK .* exp(1j*2*pi*fK*tau);

    % (3) choose DC phase = 0 or pi so H(0) is real
    phi0_est = angle(mean(H_de_K));
    phi0 = 0;
    if abs(wrapToPiLocal(phi0_est)) > (pi/2)
        phi0 = pi;
    end

    % (4) fit smooth magnitude near DC: log|H_de| vs f^2
    x = (fK.^2);
    y = log(max(abs(H_de_K), 1e-300));
    c = polyfit(x, y, polyOrder);

    % (5) synthesize fill region
    f_fill = f_pos(idx_fill);
    mag_fill = exp(polyval(c, f_fill.^2));
    H_de_fill = mag_fill .* exp(1j*phi0);
    H_fill = H_de_fill .* exp(-1j*2*pi*f_fill*tau);

    % ---- compose H_pos ----
    H_pos = zeros(size(f_pos));
    H_pos(idx_fill) = H_fill;
    H_pos(idx_meas) = H_pos_interp(idx_meas);

    % enforce DC real for real impulse
    H_pos(1) = real(H_pos(1));

    % ---- build two-sided spectrum (conjugate symmetry) ----
    Nfft_base = 2*(Npos-1);                           % even
    H_full_base = [H_pos; conj(H_pos(end-1:-1:2))];   % DC..+.. then -.. (no duplicate DC/Nyq)

    Nfft = padFactor * Nfft_base;
    nZeros = Nfft - numel(H_full_base);

    H_full = [H_full_base(1:Npos); zeros(nZeros,1); H_full_base(Npos+1:end)];

    % ---- IFFT with continuous-time scaling ----
    dt = 1/(Nfft*df);
    t  = (0:Nfft-1).' * dt;

    h = (Nfft*df) * ifft(H_full, "symmetric");
end

function x = wrapToPiLocal(x)
    x = mod((x + pi), (2*pi)) - pi;
end
