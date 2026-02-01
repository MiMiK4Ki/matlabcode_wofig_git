%% ================= main_Droute_section_run_simple_v2.m =================
% FULL EXP only / SIMPLE (no concatenation, no multi-run loop)
% - PAM/QAM supported via param.MODTYPE
% - Complex channel allowed even if PAM (no real() forcing)
% - Recommendation A applied: wTHP run re-estimates NormRef using bk

clear; close all; clc;

set(groot, 'DefaultAxesFontSize', 28);
set(groot, 'DefaultTextFontSize', 28);
set(groot, 'DefaultLegendFontSize', 28);
set(groot, 'DefaultColorbarFontSize', 28);

%% (0) parameters
param.cutoff_coeff  = 1;
param.SNRdB         = 5;

param.first_F       = -1;
param.bmax_initial  = 2;

param.MODTYPE       = "QAM";   % "PAM" or "QAM"
param.MODNUM        = 4;       % PAM: M-PAM, QAM: M-QAM (4/16/64...)
param.TXD_N         = 3000;    % 実機制約

param.hFLpS         = 50;
param.NoSpS         = 50;      % Dルートでは後で nSamps に上書き
param.fc            = 4e9;
param.beta          = 0.35;
param.Tc_but        = 1/1e9;

win = struct('first_C',-7,'C_max',15,'first_F',param.first_F,'F_max',15);

useUnitAvgPowerTx = true;   % 実験側が unit-average-power を想定なら true

%% ======= [EXP code] Path =======
addpath(genpath('C:\Users\mittaus\Desktop\Dedar'));
addpath(genpath('C:\Users\mittaus\Desktop\Dedar\iqtools'));

%% ======= [EXP code] Connection =======
Host    = 'localhost';
LANName = 'hislip0';
visaRsrc = sprintf('TCPIP0::%s::%s::INSTR', Host, LANName);

psg_device = 0;
uxr_device = 0;
vsa_soft   = 1;
connect_2_PSG_UXA_VSA;   % <-- 実験チーム側の接続関数（vsa ができる前提）

%% ======= [EXP code] HW-safe parameters (DO NOT CHANGE lightly) =======
symR        = 100e6;
frameLength = param.TXD_N;
nSamps      = 5;
Fs_bb       = symR * nSamps;

beta        = 0.35;
filtergain  = 0.4;
nsymb       = 32;

Fs_arb1     = 8e9;
fIF         = 1100e6;
gran_8190A  = 48;
tracenum    = 3;

Tc_exp = 1/symR;
fs_exp = symR*nSamps;

%% ======= my-side paramD =======
paramD = param;
paramD.NoSpS = nSamps;
paramD.beta  = beta;
paramD.fc    = fIF;
paramD.TXD_N = frameLength;

% ---- bmax の扱い（既存 thp_prepare_tx が bmax = bmax_initial*(M/2) を使う前提のまま対応）
% QAMで本来 bmax = bmax_initial*(sqrt(M)/2) を使いたい場合は、
% paramD.bmax_initial を (1/sqrt(M)) 倍しておくと、既存式で等価にできる。
modtype = upper(string(paramD.MODTYPE));
if modtype == "QAM"
    L = sqrt(double(paramD.MODNUM));
    if abs(L - round(L)) > 1e-12
        error('MODNUM must be square for QAM (4/16/64/...).');
    end
    paramD.bmax_initial = paramD.bmax_initial / L;  % ★重要（既存式のまま QAM bmax を合わせる）
end

%% ======= my-side symbols (woTHP) =======
x_train = generate_symbols(paramD, frameLength);
x_train = x_train(:).';  % row

if useUnitAvgPowerTx
    % 系列の平均電力で正規化（最小変更で安定）
    txScale = 1/sqrt(mean(abs(x_train).^2));
    x_train_tx = x_train * txScale;

    % bmaxも同じ比率で縮める（必須）
    paramD.bmax_initial = paramD.bmax_initial * txScale;
else
    x_train_tx = x_train;
end

%% -------------------------------------------------------------
%% (D1) woTHP capture
%% -------------------------------------------------------------
Swo = expteam_capture_timeseries_from_symbols_edit( ...
    x_train_tx(:), symR, nSamps, beta, filtergain, nsymb, Fs_arb1, fIF, ...
    gran_8190A, visaRsrc, tracenum, vsa);

% measured waveforms
y_rx_wo  = Swo.xout1(:).';
x_ref_wo = Swo.xin1(:).';

% ---- sync by xcorr (lag + phase) ----
c_idx_ref = (nsymb/2)*nSamps + 1;

[xc, lags] = xcorr(y_rx_wo, x_ref_wo);
[~, imx] = max(abs(xc));
lag_wo = lags(imx);
phi_wo = angle(xc(imx));

c_index_wo = c_idx_ref + lag_wo;
y_rx_wo_corr = y_rx_wo * exp(-1j*phi_wo);

% clamp (avoid out-of-range in later discretize)
c_index_wo = max(1, min(c_index_wo, numel(y_rx_wo_corr)));

fprintf('[woTHP sync] lag=%d, phi=%.3f rad, c_index=%d\n', lag_wo, phi_wo, c_index_wo);

% build IO
IO_woE = struct();
IO_woE.y_wf    = y_rx_wo_corr;
IO_woE.NoSpS   = nSamps;
IO_woE.Tc      = Tc_exp;
IO_woE.fs      = fs_exp;
IO_woE.c_index = c_index_wo;

% ---- LS estimate (use the actually sent symbol sequence) ----
coeffs_wo = estimate_channel_ls(IO_woE, x_train_tx, win);

% ---- evaluate woTHP ----
[Res_woE, D_woE] = evaluate_from_io(IO_woE, coeffs_wo, paramD, x_train_tx, 'woTHP');
evmref = NaN; if isfield(Res_woE,'EVMref'), evmref = Res_woE.EVMref; end
ser = NaN; if isfield(Res_woE,'SER'), ser = Res_woE.SER; end
fprintf('[FullEXP woTHP] BER=%.3e, SER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    Res_woE.BER, ser, Res_woE.Hard_capacity, Res_woE.Soft_capacity, Res_woE.EVM, evmref);

plot_txrx_overlay(IO_woE, coeffs_wo, paramD, x_train_tx, 'woTHP', 'D FullEXP woTHP', 200);

%% -------------------------------------------------------------
%% (D2) THP generate
%% -------------------------------------------------------------
TX_D = thp_prepare_tx(paramD, coeffs_wo, x_train_tx);

%% -------------------------------------------------------------
%% (D3) wTHP capture
%% -------------------------------------------------------------
Sw = expteam_capture_timeseries_from_symbols_edit( ...
    TX_D.bk(:), symR, nSamps, beta, filtergain, nsymb, Fs_arb1, fIF, ...
    gran_8190A, visaRsrc, tracenum, vsa);

y_rx_w  = Sw.xout1(:).';
x_ref_w = Sw.xin1(:).';

% ---- sync by xcorr (lag + phase) ----
[xc, lags] = xcorr(y_rx_w, x_ref_w);
[~, imx] = max(abs(xc));
lag_w = lags(imx);
phi_w = angle(xc(imx));

c_index_w = c_idx_ref + lag_w;
y_rx_w_corr = y_rx_w * exp(-1j*phi_w);
c_index_w = max(1, min(c_index_w, numel(y_rx_w_corr)));

fprintf('[wTHP  sync] lag=%d, phi=%.3f rad, c_index=%d\n', lag_w, phi_w, c_index_w);

IO_wE = struct();
IO_wE.y_wf    = y_rx_w_corr;
IO_wE.NoSpS   = nSamps;
IO_wE.Tc      = Tc_exp;
IO_wE.fs      = fs_exp;
IO_wE.c_index = c_index_w;

% ==========================================================
% ★ 推奨案A（必須）: wTHP run の送信列(bk)で NormRef を取り直す
% ==========================================================
coeffs_w_gain = estimate_channel_ls(IO_wE, TX_D.bk, win);

coeffs_w_eval = coeffs_wo;             % THP設計に使った係数を基本に…
coeffs_w_eval.NormRef = coeffs_w_gain.NormRef;  % NormRefだけ差し替える

ratio = coeffs_w_eval.NormRef / coeffs_wo.NormRef;
fprintf('[NormRef ratio] w/wo = %.4f%+.4fj  |.|=%.4f  ang=%.1f deg\n', ...
    real(ratio), imag(ratio), abs(ratio), angle(ratio)*180/pi);

% ---- evaluate wTHP ----
[Res_wE, D_wE] = evaluate_from_io(IO_wE, coeffs_w_eval, paramD, x_train_tx, 'wTHP');
evmref = NaN; if isfield(Res_wE,'EVMref'), evmref = Res_wE.EVMref; end
ser = NaN; if isfield(Res_wE,'SER'), ser = Res_wE.SER; end
fprintf('[FullEXP wTHP ] BER=%.3e, SER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    Res_wE.BER, ser, Res_wE.Hard_capacity, Res_wE.Soft_capacity, Res_wE.EVM, evmref);

plot_freqresp_eqchange(D_woE, D_wE, win, paramD, "D FullEXP", coeffs_wo);
plot_txrx_overlay(IO_wE, coeffs_w_eval, paramD, x_train_tx, 'wTHP', 'D FullEXP wTHP', 200, TX_D.dk);

%% ================= end =================
