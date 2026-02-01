%% ================= main_Droute_section_run_loopmean_v2.m =================
% FULL EXP only / loop & mean (NO concatenation)
% - PAM/QAM supported via param.MODTYPE
% - Complex channel allowed even if PAM
% - Recommendation A applied every run: re-estimate NormRef for wTHP using bk
% - Repeat Nrep times and print mean only

clear; close all; clc;

set(groot, 'DefaultAxesFontSize', 28);
set(groot, 'DefaultTextFontSize', 28);
set(groot, 'DefaultLegendFontSize', 28);
set(groot, 'DefaultColorbarFontSize', 28);

%% (0) parameters (edit here)
param.cutoff_coeff  = 1;
param.SNRdB         = 5;

param.first_F       = -1;
param.bmax_initial  = 2;

param.MODTYPE       = "QAM";   % "PAM" or "QAM"
param.MODNUM        = 4;       % PAM: M-PAM, QAM: M-QAM (4/16/64/...)
param.TXD_N         = 3000;    % 実機制約

param.hFLpS         = 50;
param.NoSpS         = 50;      % Dルートでは後で nSamps に上書き
param.fc            = 4e9;
param.beta          = 0.35;
param.Tc_but        = 1/1e9;

win = struct('first_C',-7,'C_max',15,'first_F',param.first_F,'F_max',15);
%% plot setting
outDir = make_run_folder(pwd, paramD, win, "Dloop", struct('symR',symR,'nSamps',nSamps,'beta',beta,'nsymb',nsymb,'filtergain',filtergain,'fIF',fIF));

saveOpt = struct('enable',true,'folder',outDir,'prefix',"Dloop",'formats',{{'png','fig'}},'dpi',200,'plotChannelCoeffTaps',true);



%% loop settings
Nrep = 10;                         % ★ここを変更
doRegenerateSymbolsEachRun = true; % true: 毎回新しいx_train, false: 固定
useUnitAvgPowerTx = false;          % true: 理想星座の平均電力で unit-average-power 化
printEachRun = true;               % 途中経過を表示
pause_between_runs_s = 0.5;        % 実機が不安定なら増やす

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
connect_2_PSG_UXA_VSA;   % vsa ができる前提

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

% ---- QAMの bmax を「既存 thp_prepare_tx の式のまま」合わせる小技 ----
% thp_prepare_tx は bmax = bmax_initial*(M/2) を使う前提。
% 望み: bmax = (user_bmax_initial)*(sqrt(M)/2)
% => paramD.bmax_initial = user_bmax_initial / sqrt(M)
% modtype = upper(string(paramD.MODTYPE));
% if modtype == "QAM"
%     L = sqrt(double(paramD.MODNUM));
%     if abs(L - round(L)) > 1e-12
%         error('MODNUM must be square for QAM (4/16/64/...).');
%     end
%     paramD.bmax_initial = paramD.bmax_initial / L;
% end

%% ---- unit-average-power scaling (theoretical, constant) ----
txScale = 1.0;
paramRun_base = paramD;

if useUnitAvgPowerTx
    M = double(paramD.MODNUM);
    if modtype == "QAM"
        L = sqrt(M);
        % square QAM (levels = -(L-1):2:(L-1)), average power = 2*(L^2-1)/3
        Pavg = 2*(L^2 - 1)/3;
    else
        % M-PAM (levels = -(M-1):2:(M-1)), average power = (M^2-1)/3
        Pavg = (M^2 - 1)/3;
    end
    txScale = 1/sqrt(Pavg);
    paramRun_base.bmax_initial = paramD.bmax_initial * txScale;
end

%% ---- allocate result arrays ----
BER_wo = nan(1,Nrep); SER_wo = nan(1,Nrep);
HC_wo  = nan(1,Nrep); SC_wo  = nan(1,Nrep);
EVM_wo = nan(1,Nrep); EVMr_wo= nan(1,Nrep);

BER_w  = nan(1,Nrep); SER_w  = nan(1,Nrep);
HC_w   = nan(1,Nrep); SC_w   = nan(1,Nrep);
EVM_w  = nan(1,Nrep); EVMr_w = nan(1,Nrep);

ratio_normref = nan(1,Nrep);   % complex ratio (w/wo)

%% ---- optional: fixed symbols if desired ----
x_train_fixed = [];
if ~doRegenerateSymbolsEachRun
    x_train_fixed = generate_symbols(paramD, frameLength);
    x_train_fixed = x_train_fixed(:).';
end

%% ================= LOOP =================
for r = 1:Nrep
    if printEachRun
        fprintf('\n================ RUN %d / %d ================\n', r, Nrep);
    end

    %% (1) symbol sequence for this run
    if doRegenerateSymbolsEachRun
        x_train = generate_symbols(paramD, frameLength);
        x_train = x_train(:).';
    else
        x_train = x_train_fixed;
    end

    x_train_tx = x_train * txScale;      % unit-average-power (optional)
    paramRun   = paramRun_base;          % includes bmax scaling if enabled

    %% -------------------------------------------------------------
    %% (D1) woTHP capture
    %% -------------------------------------------------------------
    Swo = expteam_capture_timeseries_from_symbols_edit( ...
        x_train_tx(:), symR, nSamps, beta, filtergain, nsymb, Fs_arb1, fIF, ...
        gran_8190A, visaRsrc, tracenum, vsa);

    y_rx = Swo.xout1(:).';
    x_ref = Swo.xin1(:).';

    % ---- (optional) coarse align by finddelay (keep your original style) ----
    c_idx_ref = (nsymb/2)*nSamps + 1;
    a1 = finddelay(x_ref, y_rx);
    if a1 >= 0
        y_rx = y_rx((a1+1):end);
    else
        x_ref = x_ref((-a1+1):end);
        c_idx_ref = c_idx_ref + a1; % a1 is negative -> shift ref index left
    end
    Lmin = min(numel(y_rx), numel(x_ref));
    y_rx = y_rx(1:Lmin);
    x_ref = x_ref(1:Lmin);

    % ---- xcorr sync (lag + phase) ----
    [xc, lags] = xcorr(y_rx, x_ref);
    [~, imx] = max(abs(xc));
    lag = lags(imx);
    phi = angle(xc(imx));

    c_index = c_idx_ref + lag;
    y_corr  = y_rx * exp(-1j*phi);

    c_index = max(1, min(c_index, numel(y_corr)));

    if printEachRun
        fprintf('[woTHP sync] a1=%d, lag=%d, phi=%.3f rad, c_index=%d\n', a1, lag, phi, c_index);
    end

    IO_wo = struct();
    IO_wo.y_wf    = y_corr;
    IO_wo.NoSpS   = nSamps;
    IO_wo.Tc      = Tc_exp;
    IO_wo.fs      = fs_exp;
    IO_wo.c_index = c_index;

    % LS estimate (use actually sent symbol sequence)
    coeffs_wo = estimate_channel_ls(IO_wo, x_train_tx, win);

    % evaluate woTHP
    Res_wo = evaluate_from_io(IO_wo, coeffs_wo, paramRun, x_train_tx, 'woTHP');

    BER_wo(r) = Res_wo.BER;
    if isfield(Res_wo,'SER'), SER_wo(r) = Res_wo.SER; end
    HC_wo(r)  = Res_wo.Hard_capacity;
    SC_wo(r)  = Res_wo.Soft_capacity;
    EVM_wo(r) = Res_wo.EVM;
    if isfield(Res_wo,'EVMref'), EVMr_wo(r) = Res_wo.EVMref; end

    if printEachRun
        fprintf('[RUN %d wo] BER=%.3e, SER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%\n', ...
            r, BER_wo(r), SER_wo(r), HC_wo(r), SC_wo(r), EVM_wo(r));
    end

    %% -------------------------------------------------------------
    %% (D2) THP generate
    %% -------------------------------------------------------------
    TX = thp_prepare_tx(paramRun, coeffs_wo, x_train_tx);

    %% -------------------------------------------------------------
    %% (D3) wTHP capture
    %% -------------------------------------------------------------
    Sw = expteam_capture_timeseries_from_symbols_edit( ...
        TX.bk(:), symR, nSamps, beta, filtergain, nsymb, Fs_arb1, fIF, ...
        gran_8190A, visaRsrc, tracenum, vsa);

    y_rx = Sw.xout1(:).';
    x_ref = Sw.xin1(:).';

    % coarse align by finddelay
    c_idx_ref = (nsymb/2)*nSamps + 1;
    a1 = finddelay(x_ref, y_rx);
    if a1 >= 0
        y_rx = y_rx((a1+1):end);
    else
        x_ref = x_ref((-a1+1):end);
        c_idx_ref = c_idx_ref + a1;
    end
    Lmin = min(numel(y_rx), numel(x_ref));
    y_rx = y_rx(1:Lmin);
    x_ref = x_ref(1:Lmin);

    % xcorr sync
    [xc, lags] = xcorr(y_rx, x_ref);
    [~, imx] = max(abs(xc));
    lag = lags(imx);
    phi = angle(xc(imx));

    c_index = c_idx_ref + lag;
    y_corr  = y_rx * exp(-1j*phi);
    c_index = max(1, min(c_index, numel(y_corr)));

    if printEachRun
        fprintf('[wTHP  sync] a1=%d, lag=%d, phi=%.3f rad, c_index=%d\n', a1, lag, phi, c_index);
    end

    IO_w = struct();
    IO_w.y_wf    = y_corr;
    IO_w.NoSpS   = nSamps;
    IO_w.Tc      = Tc_exp;
    IO_w.fs      = fs_exp;
    IO_w.c_index = c_index;

    % ==========================================================
    % 推奨案A（毎回適用）: wTHP run の送信列(bk)で NormRef を取り直す
    % ==========================================================
    coeffs_w_gain = estimate_channel_ls(IO_w, TX.bk, win);

    coeffs_w_eval = coeffs_wo;
    coeffs_w_eval.NormRef = coeffs_w_gain.NormRef;

    rr = coeffs_w_eval.NormRef / coeffs_wo.NormRef;
    ratio_normref(r) = rr;

    if printEachRun
        fprintf('[NormRef ratio] w/wo = %.4f%+.4fj  |.|=%.4f  ang=%.1f deg\n', ...
            real(rr), imag(rr), abs(rr), angle(rr)*180/pi);
    end

    % evaluate wTHP（ref は元の x_train_tx）
    Res_w = evaluate_from_io(IO_w, coeffs_w_eval, paramRun, x_train_tx, 'wTHP');

    BER_w(r) = Res_w.BER;
    if isfield(Res_w,'SER'), SER_w(r) = Res_w.SER; end
    HC_w(r)  = Res_w.Hard_capacity;
    SC_w(r)  = Res_w.Soft_capacity;
    EVM_w(r) = Res_w.EVM;
    if isfield(Res_w,'EVMref'), EVMr_w(r) = Res_w.EVMref; end

    if printEachRun
        fprintf('[RUN %d w ] BER=%.3e, SER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%\n', ...
            r, BER_w(r), SER_w(r), HC_w(r), SC_w(r), EVM_w(r));
    end

    pause(pause_between_runs_s);
end

%% ================= SUMMARY (MEAN ONLY) =================
fprintf('\n================ MEAN over %d runs ================\n', Nrep);

fprintf('[woTHP mean] BER=%.3e, SER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    mean(BER_wo,'omitnan'), mean(SER_wo,'omitnan'), mean(HC_wo,'omitnan'), mean(SC_wo,'omitnan'), ...
    mean(EVM_wo,'omitnan'), mean(EVMr_wo,'omitnan'));

fprintf('[wTHP  mean] BER=%.3e, SER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    mean(BER_w,'omitnan'), mean(SER_w,'omitnan'), mean(HC_w,'omitnan'), mean(SC_w,'omitnan'), ...
    mean(EVM_w,'omitnan'), mean(EVMr_w,'omitnan'));

fprintf('\n[NormRef ratio mean] mean(|w/wo|)=%.4f, mean(angle)=%.2f deg\n', ...
    mean(abs(ratio_normref),'omitnan'), mean(angle(ratio_normref)*180/pi,'omitnan'));


save(fullfile(outDir,'workspace_all.mat'),'-v7.3'); 
%% ================= end =================
