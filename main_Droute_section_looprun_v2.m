%% ================= main_Droute_section_run_v2.m =================
% FULL EXP only (D route) - PAM/QAM + complex channel + repeated runs
%
% 目的:
%   - woTHP: 実機 Tx/Rx -> 同期 -> LS推定 -> 評価
%   - wTHP : THP生成 -> 実機 Tx/Rx -> 同期 -> (推奨案A) NormRef取り直し -> 評価
%   - TXD_Nが小さい場合は、N回繰り返して統計精度を上げる
%
% 重要:
%   - PAM送信でもチャネルが複素になる可能性を許容（real() 強制なし）
%   - QAMも対応（param.MODTYPE を使用）
%   - wTHPでは runごとの未知ゲイン/位相（ARB正規化変動等）を NormRef再推定で吸収

clear; close all; clc;

set(groot, 'DefaultAxesFontSize', 28);
set(groot, 'DefaultTextFontSize', 28);
set(groot, 'DefaultLegendFontSize', 28);
set(groot, 'DefaultColorbarFontSize', 28);

%% =========================
%% (0) User parameters
%% =========================
param = struct();
param.cutoff_coeff  = 1;
param.SNRdB         = 5;

param.first_F       = -1;          % ★あなたの既存値を踏襲
param.bmax_initial  = 2;

param.MODTYPE       = "QAM";       % ★ "PAM" or "QAM"
param.MODNUM        = 4;           % ★ PAMなら M-PAM の M, QAMなら M-QAM の M(=4,16,64...)
param.TXD_N         = 3000;        % ★ 実機制約

param.hFLpS         = 50;
param.NoSpS         = 50;          % ※Dルートでは後で nSamps に上書き
param.fc            = 4e9;
param.beta          = 0.35;
param.Tc_but        = 1/1e9;
param.calc_N        = 1000;

win = struct('first_C',-7,'C_max',15,'first_F',param.first_F,'F_max',15);

%% ===== Repetition control =====
Nrep = 10;                              % ★ 例えば 10回連続測定
doRegenerateSymbolsEachRep = true;      % true: 毎回ランダムシンボル, false: 同じシンボル
useUnitAvgPowerTx = true;               % true: 送信シンボルを平均電力1にスケール（推奨）
pause_between_runs_s = 0.5;             % 実機が不安定なら増やす

%% =========================
%% EXP paths (environment)
%% =========================
addpath(genpath('C:\Users\mittaus\Desktop\Dedar'));
addpath(genpath('C:\Users\mittaus\Desktop\Dedar\iqtools'));

%% =========================
%% EXP connection (team code)
%% =========================
Host    = 'localhost';
LANName = 'hislip0';
visaRsrc = sprintf('TCPIP0::%s::%s::INSTR', Host, LANName);

psg_device = 0;
uxr_device = 0;
vsa_soft   = 1;
connect_2_PSG_UXA_VSA;   % <-- 実験チーム側の接続関数（vsa ハンドルができる前提）

%% =========================
%% HW-safe parameters (DO NOT CHANGE lightly)
%% =========================
symR        = 100e6;                 % Symbol rate [symbols/s]
frameLength = param.TXD_N;           % Symbols per frame
nSamps      = 5;                     % Samples per symbol (baseband)
Fs_bb       = symR * nSamps;
beta        = 0.35;                  % RRC rolloff
filtergain  = 0.4;                   % RRC gain
nsymb       = 32;                    % RRC span [symbols]
Fs_arb1     = 8e9;                   % DAC sample rate
fIF         = 1100e6;                % numerical upconversion center
gran_8190A  = 48;                    % AWG granularity
tracenum    = 3;                     % VSA trace index

Tc_exp = 1/symR;
fs_exp = Fs_bb;

%% =========================
%% Map to my-side paramD (this D route only)
%% =========================
paramD = param;
paramD.NoSpS = nSamps;
paramD.beta  = beta;
paramD.fc    = fIF;
paramD.TXD_N = frameLength;

%% =========================
%% Storage for per-run + concatenated evaluation
%% =========================
Res_wo_list = repmat(struct(), 1, Nrep);
Res_w_list  = repmat(struct(), 1, Nrep);
gain_ratio_list = nan(1, Nrep);

tx_cat_wo = []; rx_cat_wo = [];
tx_cat_w  = []; rx_cat_w  = [];

% optional: fix seed
% rng(1);

%% =========================
%% (D) FULL EXP loop
%% =========================
x_train_raw_cache = [];

for r = 1:Nrep
    fprintf('\n==================== RUN %d / %d ====================\n', r, Nrep);

    %% ---- Symbols (woTHP reference) ----
    if r == 1 || doRegenerateSymbolsEachRep
        x_train_raw_cache = generate_symbols(paramD, frameLength);   % PAM/QAM自動
        x_train_raw_cache = x_train_raw_cache(:).';                  % row
    end
    x_train_raw = x_train_raw_cache;

    % (optional) normalize constellation power to 1
    symScale = 1.0;
    paramRun = paramD;
    x_train_tx = x_train_raw;

    if useUnitAvgPowerTx
        Pconst = constellation_avg_power(paramRun.MODTYPE, paramRun.MODNUM);
        symScale = 1/sqrt(Pconst);
        x_train_tx = x_train_raw * symScale;

        % THP/評価で使う bmax も同じスケールに合わせる（重要）
        paramRun.bmax_initial = paramD.bmax_initial * symScale;
    end

    %% -------------------------------------------------------------
    %% (D1) woTHP capture
    %% -------------------------------------------------------------
    Swo = expteam_capture_timeseries_from_symbols_edit( ...
        x_train_tx(:), symR, nSamps, beta, filtergain, nsymb, Fs_arb1, fIF, ...
        gran_8190A, visaRsrc, tracenum, vsa);

    % build IO (xcorr sync + phase)
    [IO_woE, sync_wo] = build_IO_from_xcorr( ...
        Swo.xout1, Swo.xin1, symR, nSamps, nsymb);

    fprintf('[woTHP sync] lag=%d, phi=%.3f rad, c_index=%d\n', ...
        sync_wo.lag, sync_wo.phi, sync_wo.c_index);

    % LS estimate using the actually sent symbols (x_train_tx)
    coeffs_wo = estimate_channel_ls(IO_woE, x_train_tx, win);

    % evaluate woTHP (ref_sym must match what was actually sent)
    [Res_woE, D_woE] = evaluate_from_io(IO_woE, coeffs_wo, paramRun, x_train_tx, 'woTHP');
    Res_wo_list(r) = Res_woE;

    % concatenate aligned sequences for aggregate evaluation
    if isfield(D_woE,'tx_al') && isfield(D_woE,'rx_soft_al')
        tx_cat_wo = [tx_cat_wo, D_woE.tx_al];
        rx_cat_wo = [rx_cat_wo, D_woE.rx_soft_al];
    end

    % print per-run
    ser_wo = NaN; if isfield(Res_woE,'SER'), ser_wo = Res_woE.SER; end
    evmref = NaN; if isfield(Res_woE,'EVMref'), evmref = Res_woE.EVMref; end
    fprintf('[RUN %d woTHP] BER=%.3e, SER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
        r, Res_woE.BER, ser_wo, Res_woE.Hard_capacity, Res_woE.Soft_capacity, Res_woE.EVM, evmref);

    %% -------------------------------------------------------------
    %% (D2) THP generate (based on woTHP estimated coeffs)
    %% -------------------------------------------------------------
    TX_D = thp_prepare_tx(paramRun, coeffs_wo, x_train_tx);

    %% -------------------------------------------------------------
    %% (D3) wTHP capture
    %% -------------------------------------------------------------
    Sw = expteam_capture_timeseries_from_symbols_edit( ...
        TX_D.bk(:), symR, nSamps, beta, filtergain, nsymb, Fs_arb1, fIF, ...
        gran_8190A, visaRsrc, tracenum, vsa);

    [IO_wE, sync_w] = build_IO_from_xcorr( ...
        Sw.xout1, Sw.xin1, symR, nSamps, nsymb);

    fprintf('[wTHP  sync] lag=%d, phi=%.3f rad, c_index=%d\n', ...
        sync_w.lag, sync_w.phi, sync_w.c_index);

    % ===== 推奨案A: wTHP run の送信列(bk)で NormRef だけ取り直し =====
    [coeffs_w_eval, infoNR] = update_normref_only(coeffs_wo, IO_wE, TX_D.bk, win, "wTHP NormRef update");
    gain_ratio_list(r) = infoNR.gain_ratio;

    fprintf('[RUN %d] NormRef ratio (w/wo) = %.4f%+.4fj  |.|=%.4f  ang=%.1f deg\n', ...
        r, real(infoNR.gain_ratio), imag(infoNR.gain_ratio), abs(infoNR.gain_ratio), angle(infoNR.gain_ratio)*180/pi);

    % evaluate wTHP (ref_sym is still original x_train_tx)
    [Res_wE, D_wE] = evaluate_from_io(IO_wE, coeffs_w_eval, paramRun, x_train_tx, 'wTHP');
    Res_w_list(r) = Res_wE;

    if isfield(D_wE,'tx_al') && isfield(D_wE,'rx_soft_al')
        tx_cat_w = [tx_cat_w, D_wE.tx_al];
        rx_cat_w = [rx_cat_w, D_wE.rx_soft_al];
    end

    ser_w = NaN; if isfield(Res_wE,'SER'), ser_w = Res_wE.SER; end
    evmref = NaN; if isfield(Res_wE,'EVMref'), evmref = Res_wE.EVMref; end
    fprintf('[RUN %d wTHP ] BER=%.3e, SER=%.3e, HardCap=%.3f, SoftCap=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
        r, Res_wE.BER, ser_w, Res_wE.Hard_capacity, Res_wE.Soft_capacity, Res_wE.EVM, evmref);

    % optional plots (heavy)
    % plot_txrx_overlay(IO_woE, coeffs_wo, paramRun, x_train_tx, 'woTHP', sprintf('RUN %d woTHP',r), 200);
    % plot_txrx_overlay(IO_wE,  coeffs_w_eval, paramRun, x_train_tx, 'wTHP',  sprintf('RUN %d wTHP',r), 200, TX_D.dk);

    pause(pause_between_runs_s);
end

%% =========================
%% Aggregate evaluation (concatenated)
%% =========================
fprintf('\n==================== AGGREGATED EVALUATION ====================\n');
Agg_wo = aggregate_from_aligned(tx_cat_wo, rx_cat_wo, paramRun, symScale);
Agg_w  = aggregate_from_aligned(tx_cat_w,  rx_cat_w,  paramRun, symScale);

fprintf('[AGG woTHP] Nsym=%d, BER=%.3e, SER=%.3e, HardCap=%.3f, SoftCap(MI-hist)=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    Agg_wo.Nsym, Agg_wo.BER, Agg_wo.SER, Agg_wo.Hard_capacity, Agg_wo.Soft_capacity, Agg_wo.EVM, Agg_wo.EVMref);

fprintf('[AGG wTHP ] Nsym=%d, BER=%.3e, SER=%.3e, HardCap=%.3f, SoftCap(MI-hist)=%.3f, EVM=%.3f%%, EVMref=%.3f%%\n', ...
    Agg_w.Nsym, Agg_w.BER, Agg_w.SER, Agg_w.Hard_capacity, Agg_w.Soft_capacity, Agg_w.EVM, Agg_w.EVMref);

fprintf('\nNormRef ratio stats over runs:\n');
fprintf('  mean(|ratio|)=%.4f, std(|ratio|)=%.4f\n', mean(abs(gain_ratio_list),'omitnan'), std(abs(gain_ratio_list),'omitnan'));
fprintf('  mean(angle)=%.2f deg, std(angle)=%.2f deg\n', mean(angle(gain_ratio_list)*180/pi,'omitnan'), std(angle(gain_ratio_list)*180/pi,'omitnan'));

%% ================= end main =================


%% ==============================================================
%% Local helper functions (keep this script copy-pasteable)
%% ==============================================================

function P = constellation_avg_power(modtype, M)
% ideal constellation average power (before scaling)
    modtype = upper(string(modtype));

    if modtype == "QAM"
        L = sqrt(double(M));
        if abs(L - round(L)) > 1e-12
            error('constellation_avg_power: QAM requires square M (4,16,64,...)');
        end
        L = round(L);
        levels = double(-(L-1):2:(L-1));
        [I,Q] = meshgrid(levels, levels);
        const = I(:) + 1j*Q(:);
        P = mean(abs(const).^2);
    else
        % PAM
        levels = double(-(M-1):2:(M-1));
        P = mean(levels.^2);
    end
end

function [IO, sync] = build_IO_from_xcorr(y_rx, x_ref, symR, nSamps, nsymb)
% xcorr-based timing + phase alignment (complex OK)

    y_rx  = y_rx(:).';
    x_ref = x_ref(:).';

    Tc = 1/symR;
    fs = symR*nSamps;

    c_idx_ref = (nsymb/2)*nSamps + 1;   % group delay center (same idea as your B-route)

    [xc, lags] = xcorr(y_rx, x_ref);
    [~, imx] = max(abs(xc));

    lag = lags(imx);
    phi = angle(xc(imx));

    c_index = c_idx_ref + lag;

    % phase compensation
    y_corr = y_rx * exp(-1j*phi);

    % safety clamp
    if c_index < 1 || c_index > numel(y_corr)
        warning('build_IO_from_xcorr: c_index=%d out of range [1..%d]. Clamping.', c_index, numel(y_corr));
        c_index = min(max(c_index, 1), numel(y_corr));
    end

    IO = struct();
    IO.y_wf    = y_corr;
    IO.NoSpS   = nSamps;
    IO.Tc      = Tc;
    IO.fs      = fs;
    IO.c_index = c_index;

    sync = struct('lag',lag,'phi',phi,'c_idx_ref',c_idx_ref,'c_index',c_index);
end

function [coeffs_eval, info] = update_normref_only(coeffs_base, IO_run, tx_sym_run, win, tag)
% 推奨案A: wTHP run の送信列でLSし、NormRefのみ差し替える
    if nargin < 5, tag = ""; end

    coeffs_run = estimate_channel_ls(IO_run, tx_sym_run, win);

    coeffs_eval = coeffs_base;
    coeffs_eval.NormRef = coeffs_run.NormRef;

    info = struct();
    info.tag = string(tag);
    info.NormRef_base = coeffs_base.NormRef;
    info.NormRef_run  = coeffs_run.NormRef;

    if abs(coeffs_base.NormRef) > eps
        info.gain_ratio = coeffs_run.NormRef / coeffs_base.NormRef;
    else
        info.gain_ratio = NaN;
    end
end

function Agg = aggregate_from_aligned(tx_al, rx_al, param, symScale)
% tx_al, rx_al: symbol-rate aligned sequences (concatenated over runs)
    tx_al = tx_al(:).';
    rx_al = rx_al(:).';

    N = min(numel(tx_al), numel(rx_al));
    tx = tx_al(1:N);
    rx = rx_al(1:N);

    M = param.MODNUM;
    modtype = upper(string(param.MODTYPE));

    % ---- EVM on aligned equalized symbols ----
    EVM = NaN; EVMref = NaN;
    if exist('CalcEVMv2_silent','file') == 2
        [~,~,EVM,EVMref] = CalcEVMv2_silent(rx(:), tx(:), M);
    else
        e = rx(:) - tx(:);
        EVM = 100*sqrt(mean(abs(e).^2) / mean(abs(tx(:)).^2));
        EVMref = NaN;
    end

    % ---- Hard decisions -> SER/BER/HardCap ----
    if modtype == "QAM"
        TX_sym = qam_to_symbols(tx, M, tx);
        RX_sym = qam_to_symbols(rx, M, tx);

        TX_bits = qam_to_bits(TX_sym, M);
        RX_bits = qam_to_bits(RX_sym, M);
    else
        % PAM: rx に小さい虚部が残る場合があるので real() に寄せる
        [TX_sym, ~] = pam_to_symbols(tx, M, symScale);
        [RX_sym, ~] = pam_to_symbols(real(rx), M, symScale);

        TX_bits = pam_to_bits(tx, M, symScale);          % mappingはpam_to_bitsのデフォルト
        RX_bits = pam_to_bits(real(rx), M, symScale);
    end

    SER = mean(TX_sym ~= RX_sym);

    Lb = min(numel(TX_bits), numel(RX_bits));
    if Lb == 0
        BER = NaN;
    else
        BER = sum(xor(TX_bits(1:Lb), RX_bits(1:Lb))) / Lb;
    end

    [r_ik, HardCap] = compute_transition_and_hardcapacity(TX_sym, RX_sym, M);

    % ---- SoftCap: MI-hist estimate directly on aligned pairs ----
    SoftCap = mi_hist_aligned(tx, rx);

    Agg = struct();
    Agg.Nsym = N;
    Agg.Nb   = Lb;
    Agg.BER  = BER;
    Agg.SER  = SER;
    Agg.Hard_capacity = HardCap;
    Agg.Soft_capacity = SoftCap;
    Agg.EVM  = EVM;
    Agg.EVMref = EVMref;
    Agg.r_ik = r_ik;
end

function C = mi_hist_aligned(tx, rx)
% MI estimation by histogram (aligned pairs) - real->1D, complex->2D
% C in bits/symbol.

    tx = tx(:).';
    rx = rx(:).';

    % discrete alphabet
    symU = unique(tx);
    q = numel(symU);
    if q < 2
        C = NaN; return;
    end

    if isreal(tx) && isreal(rx)
        % 1D
        BinEdges = linspace(min(rx), max(rx), 150);
        BinWidth = BinEdges(2) - BinEdges(1);

        P = zeros(q, numel(BinEdges)-1);
        for k = 1:q
            data = rx(tx == symU(k));
            P(k,:) = histcounts(data, BinEdges, 'Normalization','pdf');
        end
        isum = sum(P, 1);

        ksum = 0;
        for k = 1:q
            pk = P(k,:);
            mask = (pk > 0) & (isum > 0);
            ksum = ksum + sum(pk(mask) .* log2(isum(mask)./pk(mask))) * BinWidth;
        end
        C = log2(q) - (1/q) * ksum;

    else
        % 2D
        Nb = 80;
        xr = real(rx); xi = imag(rx);

        BinEdgesRe = linspace(min(xr), max(xr), Nb+1);
        BinEdgesIm = linspace(min(xi), max(xi), Nb+1);
        dRe = BinEdgesRe(2) - BinEdgesRe(1);
        dIm = BinEdgesIm(2) - BinEdgesIm(1);
        BinArea = dRe*dIm;

        P = zeros(Nb, Nb, q);
        for k = 1:q
            data = rx(tx == symU(k));
            P(:,:,k) = histcounts2(real(data), imag(data), BinEdgesRe, BinEdgesIm, 'Normalization','pdf');
        end
        isum = sum(P, 3);

        ksum = 0;
        for k = 1:q
            pk = P(:,:,k);
            mask = (pk > 0) & (isum > 0);
            ksum = ksum + sum(pk(mask) .* log2(isum(mask)./pk(mask))) * BinArea;
        end
        C = log2(q) - (1/q) * ksum;
    end
end
