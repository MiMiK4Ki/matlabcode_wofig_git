%% ================= main_3routes_section_run.m =================
% 目的:
%   (1) 理論（チャネル=理論, IO=理論）で推定→THP→評価
%   (2) 一部実験（チャネル=実験推定, IO=理論）でTHP→評価
%   (3) 全部実験（チャネル=実験推定, IO=実験）で評価  ←データが揃ったらON
%
% 依存:
%   get_channel_synthetic, io_make_theory, chan_cont_to_discrete,
%   estimate_channel_ls, thp_prepare_tx, evaluate_from_io,
%   io_calc_power, io_add_awgn, generate_symbols
%
% 追加(この回答で提示):
%   chan_coeffs_to_chanframe_interp (離散coeffs→連続ChanFrame近似)

clear; close all; clc;

set(groot, 'DefaultAxesFontSize', 28);
set(groot, 'DefaultTextFontSize', 28);
set(groot, 'DefaultLegendFontSize', 28);
set(groot, 'DefaultColorbarFontSize', 28);


%% (0) 共通パラメータ（ここだけ編集）
param.cutoff_coeff  = 1;
param.SNRdB         = 10;
param.first_F       = 0;
param.bmax_initial  = 2;
param.MODNUM        = 4;           % 8-PAM (=2M)
param.MODTYPE       = "QAM";       % "PAM" or "QAM"
param.TXD_N         = 20000;

param.hFLpS         = 50;
param.NoSpS         = 50;
param.upconvert     = 10;          %#ok<NASGU>
param.fc            = 4e9;
param.beta          = 0.35;
param.Tc_but        = 1/1e9;
param.f_butterworth_cutoff = 1/param.Tc_but/2; %#ok<NASGU>
param.calc_N        = 1000; %#ok<NASGU>

win = struct('first_C',-7,'C_max',15,'first_F',param.first_F,'F_max',15);

% rng(1);  % 再現性が欲しい場合のみ

% (1) 評価に使うシンボル列（共通）
x_train = generate_symbols(param, param.TXD_N);
Nsym_train = numel(x_train);

%% ==============================================================
%% (A) 理論ルート：チャネル=理論, IO=理論
%% ==============================================================
%% (A1) 理論チャネル生成
ChanFrame_th = get_channel_synthetic(param);
figure; plot((0:numel(ChanFrame_th.h_rx)-1)*ChanFrame_th.dt, real(ChanFrame_th.h_rx));
grid on; title('ChanFrame\_th.h\_rx (real)');
figure; plot((0:numel(ChanFrame_th.h_rx)-1)*ChanFrame_th.dt, imag(ChanFrame_th.h_rx));
grid on; title('ChanFrame\_th.h\_rx (imag)');

%% (A2) 理論IO（woTHP）
IO_th_wo = io_make_theory(param, ChanFrame_th, x_train);

% タイミング（理論なので最大値でOK）
[~, IO_th_wo.c_index] = max(abs(ChanFrame_th.h_rx));
fprintf('[Theory] c_index=%d\n', IO_th_wo.c_index);

%% (A3) 雑音付加（理論：必要なら）
[P_th_wo, infoP_th_wo] = io_calc_power(IO_th_wo, Nsym_train, param, "continuous"); %#ok<NASGU>
IO_th_wo = io_add_awgn(IO_th_wo, P_th_wo, param);

%% (A4) チャネル係数（理論ブリッジ & LS推定）
coeffs_th_bridge = chan_cont_to_discrete(param, ChanFrame_th, win);
coeffs_th_ls     = estimate_channel_ls(IO_th_wo, x_train, win);

% どちらを採用するか（通常は LS を使う）
coeffs = coeffs_th_ls;
% coeffs = coeffs_th_bridge;

% 比較（任意）
figure; hold on;
stem(win.first_C:win.C_max, real(coeffs_th_bridge.ChannelCoeff), 'DisplayName','theory bridge');
stem(win.first_C:win.C_max, real(coeffs_th_ls.ChannelCoeff),     'DisplayName','LS');
grid on; legend show; title('ChannelCoeff compare (theory) [real]');
ax=gca;
ax.FontSize=28;

figure; hold on;
stem(win.first_C:win.C_max, imag(coeffs_th_bridge.ChannelCoeff), 'DisplayName','theory bridge');
stem(win.first_C:win.C_max, imag(coeffs_th_ls.ChannelCoeff),     'DisplayName','LS');
grid on; legend show; title('ChannelCoeff compare (theory) [imag]');
ax=gca;
ax.FontSize=28;

%% (A5) 評価（woTHP）
[Res_th_wo, D_th_wo] = evaluate_from_io(IO_th_wo, coeffs, param, x_train, 'woTHP');
evmref = NaN; if isfield(Res_th_wo,'EVMref'), evmref = Res_th_wo.EVMref; end
fprintf('[Theory woTHP] BER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    Res_th_wo.BER, Res_th_wo.Hard_capacity, Res_th_wo.Soft_capacity, Res_th_wo.EVM, evmref);

plot_txrx_overlay(IO_th_wo, coeffs, param, x_train, 'woTHP', 'A5 Theory woTHP', 200);


%% (A6) THP送信列 → 理論IO（wTHP）→ 評価
TX_th = thp_prepare_tx(param, coeffs, x_train);

IO_th_w = io_make_theory(param, ChanFrame_th, TX_th.bk);
IO_th_w.c_index = IO_th_wo.c_index; % とりあえず共用

Nsym_bk = numel(TX_th.bk);
[P_th_w, ~] = io_calc_power(IO_th_w, Nsym_bk, param, "continuous");
IO_th_w = io_add_awgn(IO_th_w, P_th_w, param);

% --- (A6-ADD) 推奨案A: wTHP run固有のNormRefを推定して差し替え ---
[coeffs_th_w_eval, info_th_w_norm] = coeffs_update_normref_from_io( ...
    coeffs, IO_th_w, TX_th.bk, win, "A6 Theory wTHP");

[Res_th_w, D_th_w] = evaluate_from_io(IO_th_w, coeffs_th_w_eval, param, x_train, 'wTHP');

evmref = NaN; if isfield(Res_th_w,'EVMref'), evmref = Res_th_w.EVMref; end
fprintf('[Theory wTHP ] BER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    Res_th_w.BER, Res_th_w.Hard_capacity, Res_th_w.Soft_capacity, Res_th_w.EVM, evmref);

plot_freqresp_eqchange(D_th_wo, D_th_w, win, param, "A Theory", coeffs);

plot_txrx_overlay(IO_th_w, coeffs_th_w_eval, param, x_train, 'wTHP', 'A6 Theory wTHP', 200, TX_th.dk);



%% ==============================================================
%% (B) 実験データからチャネル推定：チャネル=実験推定
%% ==============================================================
% ここは「実験データ(.mat)の変数名」に合わせて最低限だけ修正してください。
%
% 必須:
%   S.tx_upd (送信シンボル列)
%   S.xin1   (Tx側で既知の連続参照波形)  ※無ければ RaisedCosineTransmitFilter で生成
%   S.xout1  (Rx連続波形)
%   S.symR, S.nSamps
%
% ※ここでは「推定した coeffs_exp を後段の(一部実験/全部実験)に渡す」ことだけやる。

%% (B1) 実験データ load
% measFile = 'IQWAVEFORM_64QAM.mat';
measFile = 'IQWAVEFORM_64QAM.mat';   % ←適宜変更
S = load(measFile);

tx_sym_exp = S.tx_upd(:).';   % 送信シンボル列（実験）
y_rx_exp   = S.xout1(:).';    % 受信波形（実験）
x_ref_exp  = S.xin1(:).';     % 参照（実験送信側で既知の連続波形）

symR_exp   = S.symR;
NoSpS_exp  = S.nSamps;

Tc_exp = 1/symR_exp;
fs_exp = symR_exp*NoSpS_exp;

% RC Tx filter パラメータ（ファイルに無ければ固定値）
beta = 0.35;
nsymb = 80;
% filtergain は実験コードに合わせる（ファイルに無ければ固定値）
filtergain = 0.3578;

%% (B2) 参照側の記号中心（群遅延）c_idx_ref
c_idx_ref = (nsymb/2)*NoSpS_exp + 1;
y_rx_exp = [zeros(1, 10) y_rx_exp]
%% (B3) xcorr による lag と位相推定 → c_index を決める
[xc, lags] = xcorr(y_rx_exp, x_ref_exp);
[~, imx] = max(abs(xc));
lag = lags(imx);
phi = angle(xc(imx));

c_index_exp = c_idx_ref + lag;
fprintf('[EXP sync] lag=%d, phi=%.3f rad, c_index=%d\n', lag, phi, c_index_exp);

% 受信位相補償（参照に合わせる）
y_rx_exp_corr = y_rx_exp * exp(-1j*phi);

%% (B4) 実験IOを作る（推定用）
IO_exp = struct();
IO_exp.y_wf    = y_rx_exp_corr;
IO_exp.NoSpS   = NoSpS_exp;
IO_exp.Tc      = Tc_exp;
IO_exp.fs      = fs_exp;
IO_exp.c_index = c_index_exp;

%% (B5) 実験IOから LS で coeffs を推定
coeffs_exp = estimate_channel_ls(IO_exp, tx_sym_exp, win);

% =========================================================
% 追加: Tx RaisedCosineTransmitFilter の連続インパルス応答を重ねて確認
%   - 実験側コードと同じ comm.RaisedCosineTransmitFilter を使う
%   - 連続インパルス応答を「ピーク位置」で正規化して、離散推定 taps と比較する
%   ※注意: 推定 taps は「Txフィルタ＋チャネル＋Rx」を含む可能性があるので、
%          完全一致しなくてもOK（形の確認が主目的）
% =========================================================

rcFilt_check = comm.RaisedCosineTransmitFilter( ...
    'RolloffFactor', beta, ...
    'FilterSpanInSymbols', nsymb, ...
    'OutputSamplesPerSymbol', NoSpS_exp, ...
    'Gain', filtergain);

reset(rcFilt_check);
h_rc = step(rcFilt_check, [1; zeros(nsymb,1)]);   % 連続インパルス応答（サンプル列）

% ピークで中心合わせ
[~, c_idx_peak_rc] = max(abs(h_rc));
h_rc_n = h_rc / h_rc(c_idx_peak_rc) .* coeffs_exp.ChannelCoeff(1-win.first_C) ;              % ピークをシンボルと合わせる

% 横軸を「シンボルインデックス相当」に変換（連続）
m_frac = ((0:numel(h_rc_n)-1) - (c_idx_peak_rc-1)) / NoSpS_exp;

% 離散 taps（推定値）は m=win.first_C:win.C_max
m_int = win.first_C:win.C_max;

figure; hold on;
plot(m_frac, real(h_rc_n), 'DisplayName','|Tx RC impulse| (real)');
stem(m_int, real(coeffs_exp.ChannelCoeff), 'filled', 'DisplayName','|Estimated ChannelCoeff| (real)');
grid on;
xlabel('symbol index m'); ylabel('magnitude');
title('B5 Tx RC impulse vs LS-estimated channel [real]');
legend('show','Location','best');

figure; hold on;
plot(m_frac, imag(h_rc_n), 'DisplayName','|Tx RC impulse| (imag)');
stem(m_int, imag(coeffs_exp.ChannelCoeff), 'filled', 'DisplayName','|Estimated ChannelCoeff| (imag)');
grid on;
xlabel('symbol index m'); ylabel('magnitude');
title('B5 Tx RC impulse vs LS-estimated channel [imag]');
legend('show','Location','best');

ax=gca;
ax.FontSize=28;
%% ==============================================================
%% (C) 一部実験ルート：チャネル=実験推定, IO=理論（シミュレーション）
%% ==============================================================
% ここでは coeffs_exp を使って、連続IRを近似生成し、io_make_theory で通す。
% 近似なので「連続波形の形」は厳密一致しないが、
%   - サンプリング点の離散モデル
%   - THPの符号化/復号
% の動作確認が主目的。
%% チャネル係数の実数化（PAM用。QAMの場合は複素を保持）
% if upper(string(param.MODTYPE)) == "PAM"
%     coeffs_exp.ChannelCoeff = real(coeffs_exp.ChannelCoeff);
%     coeffs_exp.FilterCoeff = real(coeffs_exp.FilterCoeff);
%     coeffs_exp.NormRef = real(coeffs_exp.NormRef);
% end

coeffs_exp.ChannelCoeff(1:-win.first_C-1) = 0;
coeffs_exp.FilterCoeff(1:-win.first_F-1) = 0;
%% (C1) 実験推定 coeffs → 近似 ChanFrame を作る（連続IRを生成）
ChanFrame_px = chan_coeffs_to_chanframe_interp(param, coeffs_exp, win, "linear");

% この近似ChanFrameで、逆に離散タップを取り出して一致チェック（推奨）
coeffs_px_check = chan_cont_to_discrete(param, ChanFrame_px, win);
figure; hold on;
stem(win.first_C:win.C_max, real(coeffs_exp.ChannelCoeff),      'DisplayName','coeffs\_exp');
stem(win.first_C:win.C_max, real(coeffs_px_check.ChannelCoeff), 'DisplayName','re-extract');
plot(win.first_C:1/param.NoSpS:(win.C_max+5),real(ChanFrame_px.h_rx ./ coeffs_exp.NormRef),'DisplayName','Interp (real)')
grid on; legend show; title('Sanity: coeffs_exp vs re-extracted [real]');

figure; hold on;
stem(win.first_C:win.C_max, imag(coeffs_exp.ChannelCoeff),      'DisplayName','coeffs\_exp');
stem(win.first_C:win.C_max, imag(coeffs_px_check.ChannelCoeff), 'DisplayName','re-extract');
plot(win.first_C:1/param.NoSpS:(win.C_max+5),imag(ChanFrame_px.h_rx ./ coeffs_exp.NormRef),'DisplayName','Interp (imag)')
grid on; legend show; title('Sanity: coeffs_exp vs re-extracted [imag]');

ax=gca;
ax.FontSize=28;

%% (C2) 一部実験：IO(woTHP)を理論で生成 → 評価
IO_px_wo = io_make_theory(param, ChanFrame_px, x_train);

% この近似ChanFrameでは m=0 の位置を明示的に中心にする
IO_px_wo.c_index = 1 + (0 - win.first_C)*param.NoSpS;

[P_px_wo, ~] = io_calc_power(IO_px_wo, Nsym_train, param, "continuous");
IO_px_wo = io_add_awgn(IO_px_wo, P_px_wo, param);   % 加えたくなければコメントアウト

[Res_px_wo, D_px_wo] = evaluate_from_io(IO_px_wo, coeffs_exp, param, x_train, 'woTHP');

evmref = NaN; if isfield(Res_px_wo,'EVMref'), evmref = Res_px_wo.EVMref; end
fprintf('[PartialExp woTHP] BER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    Res_px_wo.BER, Res_px_wo.Hard_capacity, Res_px_wo.Soft_capacity, Res_px_wo.EVM, evmref);

plot_txrx_overlay(IO_px_wo, coeffs_exp, param, x_train, 'woTHP', 'C2 PartialExp woTHP', 200);



%% (C3) 一部実験：THP → 理論IO → 評価
TX_px = thp_prepare_tx(param, coeffs_exp, x_train);

IO_px_w = io_make_theory(param, ChanFrame_px, TX_px.bk);
IO_px_w.c_index = IO_px_wo.c_index;

Nsym_bk = numel(TX_px.bk);
[P_px_w, ~] = io_calc_power(IO_px_w, Nsym_bk, param, "continuous");
IO_px_w = io_add_awgn(IO_px_w, P_px_w, param);      % 加えたくなければコメントアウト

% --- (C3-ADD) 推奨案A: wTHP run固有のNormRefを推定して差し替え ---
[coeffs_px_w_eval, info_px_w_norm] = coeffs_update_normref_from_io( ...
    coeffs_exp, IO_px_w, TX_px.bk, win, "C3 PartialExp wTHP");

[Res_px_w, D_px_w] = evaluate_from_io(IO_px_w, coeffs_px_w_eval, param, x_train, 'wTHP');

evmref = NaN; if isfield(Res_px_w,'EVMref'), evmref = Res_px_w.EVMref; end
fprintf('[PartialExp wTHP ] BER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    Res_px_w.BER, Res_px_w.Hard_capacity, Res_px_w.Soft_capacity, Res_px_w.EVM, evmref);

plot_freqresp_eqchange(D_px_wo, D_px_w, win, param, "C PartialExp", coeffs_exp);

plot_txrx_overlay(IO_px_w, coeffs_px_w_eval, param, x_train, 'wTHP', 'C3 PartialExp wTHP', 200, TX_px.dk);


%% ==============================================================
%% (E) 一部実験ルート（S2P）：チャネル=s2p(デバイス) + sRRC, IO=理論
%% ==============================================================
s2pFile = "Filter response.s2p";   % ★適宜変更
Zs = 50; Zl = 50;

optS2P = struct();
optS2P.Kfit      = 25;      % -3dB基準/低域位相推定に使う点数（小さすぎ注意）
optS2P.doFillLF  = false;    % f< f_start の低域補完をする
optS2P.doPlot    = true;    % H(f), h(t) 可視化
optS2P.padPow2   = true;    % Nfft を 2^n に丸める

%% --- (E1) s2p + sRRC から ChanFrame を作る（A1と同じ conv*dt 規約） ---
[ChanFrame_s2p, info_s2p, param_s2p] = get_channel_s2p_device_srrc(param, win, s2pFile, Zs, Zl, optS2P);

fprintf('[S2P] Bdev(-3dB)=%.3f GHz, Bsrrc=Bdev*cutoff=%.3f GHz\n', ...
    info_s2p.Bdev_Hz/1e9, info_s2p.Bsrrc_Hz/1e9);
fprintf('[S2P] Tc=%.3f ps, dt=%.3f ps, fs=%.3f GHz\n', ...
    info_s2p.Tc*1e12, info_s2p.dt*1e12, info_s2p.fs/1e9);

%% --- (E2) IO(woTHP) ---
IO_s2p_wo = io_make_theory(param_s2p, ChanFrame_s2p, x_train);
[~, IO_s2p_wo.c_index] = max(abs(ChanFrame_s2p.h_rx));
fprintf('[S2P] c_index=%d\n', IO_s2p_wo.c_index);

[P_s2p_wo, ~] = io_calc_power(IO_s2p_wo, Nsym_train, param_s2p, "continuous");
IO_s2p_wo = io_add_awgn(IO_s2p_wo, P_s2p_wo, param_s2p);

%% --- (E3) 係数は "連続→離散ブリッジ" を主に使う（s2pは既知チャネルなので） ---
coeffs_s2p = chan_cont_to_discrete(param_s2p, ChanFrame_s2p, win);

%% （任意）LSと比較したければ

peak_shift_index = IO_s2p_wo.c_index - (2*param.hFLpS*param.NoSpS + 1);

coeffs_s2p_ls = estimate_channel_ls(IO_s2p_wo, x_train, win);
figure; hold on;
plot((-2*param.hFLpS :1/param.NoSpS: 2*param.hFLpS)  - peak_shift_index/param.NoSpS   ,real(ChanFrame_s2p.h_rx ./coeffs_s2p.NormRef), 'DisplayName','Continuous (real)')
stem(win.first_C:win.C_max, real(coeffs_s2p.ChannelCoeff), 'filled', 'DisplayName','bridge (real)');
stem(win.first_C:win.C_max, real(coeffs_s2p_ls.ChannelCoeff),     'DisplayName','LS (real)');
grid on; legend show; title('S2P route ChannelCoeff: bridge vs LS [real]');
xlim([win.first_C win.C_max])

figure; hold on;
plot((-2*param.hFLpS :1/param.NoSpS: 2*param.hFLpS)  - peak_shift_index/param.NoSpS   ,imag(ChanFrame_s2p.h_rx ./coeffs_s2p.NormRef), 'DisplayName','Continuous (imag)')
stem(win.first_C:win.C_max, imag(coeffs_s2p.ChannelCoeff), 'filled', 'DisplayName','bridge (imag)');
stem(win.first_C:win.C_max, imag(coeffs_s2p_ls.ChannelCoeff),     'DisplayName','LS (imag)');
grid on; legend show; title('S2P route ChannelCoeff: bridge vs LS [imag]');
xlim([win.first_C win.C_max])
%% --- (E4) 評価 woTHP ---
[Res_s2p_wo, D_s2p_wo] = evaluate_from_io(IO_s2p_wo, coeffs_s2p, param_s2p, x_train, 'woTHP');
evmref = NaN; if isfield(Res_s2p_wo,'EVMref'), evmref = Res_s2p_wo.EVMref; end
fprintf('[S2P woTHP] BER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    Res_s2p_wo.BER, Res_s2p_wo.Hard_capacity, Res_s2p_wo.Soft_capacity, Res_s2p_wo.EVM, evmref);

plot_txrx_overlay(IO_s2p_wo, coeffs_s2p, param_s2p, x_train, 'woTHP', 'E2 S2P woTHP', 200);

%% --- (E5) THP -> IO(wTHP) ---
TX_s2p = thp_prepare_tx(param_s2p, coeffs_s2p, x_train);

IO_s2p_w = io_make_theory(param_s2p, ChanFrame_s2p, TX_s2p.bk);
IO_s2p_w.c_index = IO_s2p_wo.c_index;

Nsym_bk_s2p = numel(TX_s2p.bk);
[P_s2p_w, ~] = io_calc_power(IO_s2p_w, Nsym_bk_s2p, param_s2p, "continuous");
IO_s2p_w = io_add_awgn(IO_s2p_w, P_s2p_w, param_s2p);

%% --- (E6) 評価 wTHP ---
% --- (E6-ADD) 推奨案A: wTHP run固有のNormRefを推定して差し替え ---
[coeffs_s2p_w_eval, info_s2p_w_norm] = coeffs_update_normref_from_io( ...
    coeffs_s2p, IO_s2p_w, TX_s2p.bk, win, "E6 S2P wTHP");

[Res_s2p_w, D_s2p_w] = evaluate_from_io(IO_s2p_w, coeffs_s2p_w_eval, param_s2p, x_train, 'wTHP');

evmref = NaN; if isfield(Res_s2p_w,'EVMref'), evmref = Res_s2p_w.EVMref; end
fprintf('[S2P wTHP ] BER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    Res_s2p_w.BER, Res_s2p_w.Hard_capacity, Res_s2p_w.Soft_capacity, Res_s2p_w.EVM, evmref);

[h_wo,h_w] = plot_freqresp_eqchange(D_s2p_wo, D_s2p_w, win, param_s2p, "E S2P", coeffs_s2p);

lin_conv_output = conv(h_w,x_train);
lin_conv_output = lin_conv_output(abs(win.first_C)+1:200 + abs(win.first_C));

plot_txrx_overlay(IO_s2p_w, coeffs_s2p_w_eval, param_s2p, x_train, 'wTHP', 'E3 S2P wTHP', 200, TX_s2p.dk, lin_conv_output);

%%
[h_wo,h_w] = plot_freqresp_eqchange(D_s2p_wo, D_s2p_w, win, param_s2p, "E S2P", coeffs_s2p);

lin_conv_output = conv(h_w,x_train);
lin_conv_output = lin_conv_output(abs(win.first_C)+1:200 + abs(win.first_C));


%%
plot_txrx_overlay(IO_s2p_w, coeffs_s2p, param_s2p, x_train, 'wTHP', 'E3 S2P wTHP', 200, TX_s2p.dk,lin_conv_output);

%% ==============================================================
%% (D) 全部実験ルート：チャネル=実験推定, IO=実験（save/load無しの一体運用）
%% ==============================================================
% 流れ:
%   1) woTHP実験測定（tx_upd = x_train） → IO作成 → coeffs_exp推定 → woTHP評価
%   2) THP生成 → wTHP実験測定（tx_upd = TX.bk） → IO作成 → wTHP評価
%
% ※ 実機が繋がっているPCでのみ動作
% ※ Dedar/iqtools のパスは環境に合わせて要修正

doRouteD_fullExperiment = true;

if doRouteD_fullExperiment

    %% (D0) 実験設定（EXPチームコード由来の変数名を維持）
    expCfg = struct();

    expCfg.psg_device = 0;
    expCfg.uxr_device = 1;
    expCfg.vsa_soft   = 1;

    expCfg.paths = { ...
        'C:\Users\mittaus\Desktop\Dedar', ...
        'C:\Users\mittaus\Desktop\iqtools'};

    expCfg.symR       = 1000e6;
    expCfg.nSamps     = 8;

    expCfg.beta       = 0.35;
    expCfg.filtergain = 0.3578;
    expCfg.nsymb      = 80;

    expCfg.Fs_ARB   = 128e9;   % ARB側のサンプルレート（実験コードの想定）
    expCfg.fc_shift = 4e9;     % iqloadfileで使ってたfrequencyShift相当

    expCfg.low_voltage_levels  = -0.83;
    expCfg.high_voltage_levels =  0.83;

    expCfg.tracenum = [3];
    expCfg.capture_pause = 5;         % 秒
    expCfg.normalize_rx_rms = true;   % 実験コードの xout1=rms正規化を踏襲

    expCfg.use_iqloadfile_tempfile = true;   % 最小変更で動きやすい（関数内でtemp保存→削除）

    %% (D0-1) 実験機器接続（EXPチームコード由来）
    EXP = exp_open(expCfg);

    try
        %% --------------------------
        %% (D1) woTHP：実験測定 → 推定 → 評価
        %% --------------------------
        tx_upd = x_train(:);   % ★ここが「実験チームのqammod部分」を置換するあなた側の本体

        % --- 実験測定（EXPチームコードを関数化したもの） ---
        Swo = exp_measure_once(EXP, tx_upd, expCfg);

        % --- (B2)-(B4) あなた側ロジックで同期して IO化 ---
        [IO_woE, lag_wo, phi_wo, c_index_wo] = exp_build_IO_from_S(Swo);

        fprintf('[FullEXP woTHP sync] lag=%d, phi=%.3f rad, c_index=%d\n', lag_wo, phi_wo, c_index_wo);

        % --- param を実験条件に合わせたコピー（安全策） ---
        param_exp = param;
        param_exp.NoSpS = Swo.nSamps;

        % --- (B5) LSでcoeffs推定（実験データ） ---
        coeffs_exp = estimate_channel_ls(IO_woE, tx_upd(:).', win);

        % Cルート等で破壊されないように退避（あなたの元コード事情）
        coeffs_exp_est = coeffs_exp;

        % --- woTHP評価（実験IOで評価） ---
        [Res_woE, D_woE] = evaluate_from_io(IO_woE, coeffs_exp_est, param_exp, x_train, 'woTHP');

        evmref = NaN; if isfield(Res_woE,'EVMref'), evmref = Res_woE.EVMref; end
        fprintf('[FullEXP woTHP] BER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
            Res_woE.BER, Res_woE.Hard_capacity, Res_woE.Soft_capacity, Res_woE.EVM, evmref);

        plot_txrx_overlay(IO_woE, coeffs_exp_est, param_exp, x_train, 'woTHP', 'D1 FullExp woTHP', 200);

        %% --------------------------
        %% (D2) THP生成（あなた側）
        %% --------------------------
        TX_D = thp_prepare_tx(param_exp, coeffs_exp_est, x_train);

        %% --------------------------
        %% (D3) wTHP：実験測定 → 評価
        %% --------------------------
        tx_upd_w = TX_D.bk(:);   % ★これが "推定→THP→実験送信" の肝

        Sw = exp_measure_once(EXP, tx_upd_w, expCfg);

        [IO_wE, lag_w, phi_w, c_index_w] = exp_build_IO_from_S(Sw);
        fprintf('[FullEXP wTHP  sync] lag=%d, phi=%.3f rad, c_index=%d\n', lag_w, phi_w, c_index_w);

        param_exp_w = param_exp;
        param_exp_w.NoSpS = Sw.nSamps;

        % --- (D3-ADD) 推奨案A: wTHP run固有のNormRefを推定して差し替え ---
        [coeffs_wE_eval, info_wE_norm] = coeffs_update_normref_from_io( ...
            coeffs_exp_est, IO_wE, tx_upd_w(:).', win, "D3 FullExp wTHP");

        [Res_wE, D_wE] = evaluate_from_io(IO_wE, coeffs_wE_eval, param_exp_w, x_train, 'wTHP');

        evmref = NaN; if isfield(Res_wE,'EVMref'), evmref = Res_wE.EVMref; end
        fprintf('[FullEXP wTHP ] BER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
            Res_wE.BER, Res_wE.Hard_capacity, Res_wE.Soft_capacity, Res_wE.EVM, evmref);

        plot_freqresp_eqchange(D_woE, D_wE, win, param_exp, "D FullExp", coeffs_exp_est);

        plot_txrx_overlay(IO_wE, coeffs_wE_eval, param_exp_w, x_train, 'wTHP', 'D2 FullExp wTHP', 200, TX_D.dk);


    catch ME
        exp_close(EXP);
        rethrow(ME);
    end

    exp_close(EXP);
end

%% ================= end =================
