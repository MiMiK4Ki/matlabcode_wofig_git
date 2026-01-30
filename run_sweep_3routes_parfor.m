%% ================= run_sweep_3routes_parfor.m =================
% 複数パラメータを parfor で回す版（プロット無し）
%
% 追加:
%   - s2p ルート (device s2p + sRRC) を追加して保存
%
% IMPORTANT(TODO):
%   いまは実験coeffsの複素回転が未対応なので、暫定的に real(coeffs_exp) で動かす場合がある。
%   後で必ず「A0による位相補償＋実数軸投影」に戻すこと。  ←忘れないこと！

clear; close all; clc;

%% -----------------------
%% (0) Sweepしたいパラメータ（ここだけ編集）
%% -----------------------
cutoff_list = [0.6 0.8 1.0 1.2 1.4 1.6];
snr_list    = [0:2:30];
firstF_list = [-2 -1 0 1 2];
bmax_list   = [2:1:15];
mod_list    = [2 4 8 16];     % (=MODNUM)

%% ★S2P: ルート設定
useS2P  = true;
s2pFile = "Filter response.s2p";   % ★適宜変更
Zs = 50; Zl = 50;

optS2P = struct();
optS2P.Kfit     = 25;
optS2P.doFillLF = false;
optS2P.doPlot   = false;      % sweepでは必ずfalse（parfor前提）
optS2P.padPow2  = true;

%% -----------------------
%% (1) 固定パラメータ（共通）
%% -----------------------
baseParam = struct();
baseParam.cutoff_coeff  = cutoff_list(1);
baseParam.SNRdB         = snr_list(1);
baseParam.first_F       = firstF_list(1);
baseParam.bmax_initial  = bmax_list(1);
baseParam.MODNUM        = mod_list(1);
baseParam.MODTYPE       = "PAM";        % "PAM" or "QAM"
baseParam.TXD_N         = 20090;

baseParam.hFLpS         = 50;
baseParam.NoSpS         = 50;
baseParam.upconvert     = 10;          %#ok<NASGU>
baseParam.fc            = 4e9;
baseParam.beta          = 0.35;

baseParam.Tc_but        = 1/1e9;
baseParam.f_butterworth_cutoff = 1/baseParam.Tc_but/2; %#ok<NASGU>
baseParam.calc_N        = 1000; %#ok<NASGU>

% winの固定部（first_Fだけケースごとに入れ替え）
winBase = struct('first_C',-3,'C_max',15,'first_F',0,'F_max',15);

%% -----------------------
%% (2) 事前計算（Broadcast化して parfor の中を軽くする）
%% -----------------------

% 2-1) x_train を MODNUMごとに作る（SNRやcutoffが違っても同じ列を使う）
x_train_by_mod = cell(1, numel(mod_list));
rng(1);
for im = 1:numel(mod_list)
    tmpParam = baseParam;
    tmpParam.MODNUM = mod_list(im);
    x_train_by_mod{im} = generate_symbols(tmpParam, baseParam.TXD_N);
end

% 2-2) 理論ChanFrameを cutoffごとに作る（cutoff以外が固定ならこれでOK）
ChanFrame_th_by_cutoff = cell(1, numel(cutoff_list));
for ic = 1:numel(cutoff_list)
    ptmp = baseParam;
    ptmp.cutoff_coeff = cutoff_list(ic);
    ChanFrame_th_by_cutoff{ic} = get_channel_synthetic(ptmp);
end

% ★S2P 2-3) s2p ChanFrame を cutoff ごとに作る（device+sRRC 合成済み）
ChanFrame_s2p_by_cutoff = cell(1, numel(cutoff_list));
param_s2p_by_cutoff     = cell(1, numel(cutoff_list));  % Tc_but 等が上書きされた param_out
info_s2p_by_cutoff      = cell(1, numel(cutoff_list));

if useS2P
    if ~exist(s2pFile,'file')
        error('S2P file not found: %s', s2pFile);
    end

    for ic = 1:numel(cutoff_list)
        ptmp = baseParam;
        ptmp.cutoff_coeff = cutoff_list(ic);

        % get_channel_s2p_device_srrc は「conv(h_dev, h_srrc)*dt」規約で ChanFrame を返す前提
        [ChanFrame_s2p_by_cutoff{ic}, info_s2p_by_cutoff{ic}, param_s2p_by_cutoff{ic}] = ...
            get_channel_s2p_device_srrc(ptmp, winBase, s2pFile, Zs, Zl, optS2P);
    end
end

%% -----------------------
%% (3) 実験チャネル推定（必要なら先に1回だけ作る）
%% -----------------------
usePartialExp = true;               % ←不要なら false
coeffs_exp_by_firstF = [];

if usePartialExp
    measFile = 'IQWAVEFORM_64QAM.mat';   % 適宜
    S = load(measFile);

    tx_sym_exp = S.tx_upd(:).';     % 実験送信シンボル（QAMで複素になり得る）
    y_rx_exp   = S.xout1(:).';
    x_ref_exp  = S.xin1(:).';

    symR_exp   = S.symR;
    NoSpS_exp  = S.nSamps;

    % 群遅延（実験側コードに合わせる）
    beta = 0.35; %#ok<NASGU>
    nsymb = 80;  %#ok<NASGU>
    c_idx_ref = (nsymb/2)*NoSpS_exp + 1;

    % xcorr で lag & 位相 -> c_index
    [xc, lags] = xcorr(y_rx_exp, x_ref_exp);
    [~, imx]   = max(abs(xc));
    lag        = lags(imx);
    phi        = angle(xc(imx));
    c_index_exp = c_idx_ref + lag;

    y_rx_exp_corr = y_rx_exp * exp(-1j*phi);

    IO_exp = struct();
    IO_exp.y_wf    = y_rx_exp_corr;
    IO_exp.NoSpS   = NoSpS_exp;
    IO_exp.Tc      = 1/symR_exp;
    IO_exp.fs      = symR_exp*NoSpS_exp;
    IO_exp.c_index = c_index_exp;

    % first_Fごとに coeffs を作る（sweepでfirst_Fを変えるならこれが簡単）
    coeffs_exp_by_firstF = cell(1, numel(firstF_list));
    for iF = 1:numel(firstF_list)
        winTmp = winBase;
        winTmp.first_F = firstF_list(iF);
        coeffs_tmp = estimate_channel_ls(IO_exp, tx_sym_exp, winTmp);

        % -------------------------
        % 暫定：複素を real に丸める（TODO: 後で必ず直す）
        % -------------------------
        if upper(string(baseParam.MODTYPE)) == "PAM"
            if any(abs(imag(coeffs_tmp.ChannelCoeff)) > 1e-12) || abs(imag(coeffs_tmp.NormRef)) > 1e-12
                warning('coeffs_exp is complex. PAM mode: using real(coeffs_exp).');
            end
            coeffs_tmp.ChannelCoeff = real(coeffs_tmp.ChannelCoeff);
            coeffs_tmp.FilterCoeff  = real(coeffs_tmp.FilterCoeff);
            coeffs_tmp.NormRef      = real(coeffs_tmp.NormRef);
        end

        coeffs_exp_by_firstF{iF} = coeffs_tmp;
    end
end

%% -----------------------
%% (4) パラメータグリッド（idx=1..meshN の整数連続に落とす）
%% -----------------------
[nC, nS, nF, nB, nM] = deal(numel(cutoff_list), numel(snr_list), numel(firstF_list), numel(bmax_list), numel(mod_list));
gridbox = [nC nS nF nB nM];

[cutI, snrI, fFI, bmaxI, modI] = ndgrid(1:nC, 1:nS, 1:nF, 1:nB, 1:nM);
meshN = numel(cutI);

cutI  = cutI(:);
snrI  = snrI(:);
fFI   = fFI(:);
bmaxI = bmaxI(:);
modI  = modI(:);

%% -----------------------
%% (5) parfor で回す（プロット禁止）
%% -----------------------
results = cell(meshN, 1);

useParfor = false;
if useParfor
    parfor idx = 1:meshN
        rng(1000 + idx);

        param = baseParam;
        param.cutoff_coeff = cutoff_list(cutI(idx));
        param.SNRdB        = snr_list(snrI(idx));
        param.first_F      = firstF_list(fFI(idx));
        param.bmax_initial = bmax_list(bmaxI(idx));
        param.MODNUM       = mod_list(modI(idx));

        win = winBase;
        win.first_F = param.first_F;

        x_train = x_train_by_mod{modI(idx)};
        ChanFrame_th = ChanFrame_th_by_cutoff{cutI(idx)};

        if usePartialExp
            coeffs_exp = coeffs_exp_by_firstF{fFI(idx)};
        else
            coeffs_exp = [];
        end

        % ★S2P: cutoffに対応する ChanFrame と param_out を渡す
        if useS2P
            ChanFrame_s2p = ChanFrame_s2p_by_cutoff{cutI(idx)};
            param_s2p_base = param_s2p_by_cutoff{cutI(idx)};
        else
            ChanFrame_s2p = [];
            param_s2p_base = [];
        end

        results{idx} = run_one_case_3routes_noplot(param, win, ChanFrame_th, x_train, coeffs_exp, ChanFrame_s2p, param_s2p_base);
    end
else
    for idx = 1:meshN
        rng(1000 + idx);

        param = baseParam;
        param.cutoff_coeff = cutoff_list(cutI(idx));
        param.SNRdB        = snr_list(snrI(idx));
        param.first_F      = firstF_list(fFI(idx));
        param.bmax_initial = bmax_list(bmaxI(idx));
        param.MODNUM       = mod_list(modI(idx));

        win = winBase;
        win.first_F = param.first_F;

        x_train = x_train_by_mod{modI(idx)};
        ChanFrame_th = ChanFrame_th_by_cutoff{cutI(idx)};

        if usePartialExp
            coeffs_exp = coeffs_exp_by_firstF{fFI(idx)};
        else
            coeffs_exp = [];
        end

        % ★S2P
        if useS2P
            ChanFrame_s2p = ChanFrame_s2p_by_cutoff{cutI(idx)};
            param_s2p_base = param_s2p_by_cutoff{cutI(idx)};
        else
            ChanFrame_s2p = [];
            param_s2p_base = [];
        end

        results{idx} = run_one_case_3routes_noplot(param, win, ChanFrame_th, x_train, coeffs_exp, ChanFrame_s2p, param_s2p_base);
    end
end

%% -----------------------
%% (6) 保存（Param + Res だけ）  ※EVM/EVMref/S2Pも追加
%% -----------------------
R = [results{:}];

res = struct();

% Theory route
res.th_BER_woTHP     = reshape([R.th_BER_woTHP],     gridbox);
res.th_BER_wTHP      = reshape([R.th_BER_wTHP],      gridbox);
res.th_HardCap_woTHP = reshape([R.th_HardCap_woTHP], gridbox);
res.th_HardCap_wTHP  = reshape([R.th_HardCap_wTHP],  gridbox);
res.th_SoftCap_woTHP = reshape([R.th_SoftCap_woTHP], gridbox);
res.th_SoftCap_wTHP  = reshape([R.th_SoftCap_wTHP],  gridbox);
res.th_EVM_woTHP     = reshape([R.th_EVM_woTHP],     gridbox);
res.th_EVM_wTHP      = reshape([R.th_EVM_wTHP],      gridbox);
res.th_EVMref_woTHP  = reshape([R.th_EVMref_woTHP],  gridbox);
res.th_EVMref_wTHP   = reshape([R.th_EVMref_wTHP],   gridbox);

% PartialExp route
res.px_BER_woTHP     = reshape([R.px_BER_woTHP],     gridbox);
res.px_BER_wTHP      = reshape([R.px_BER_wTHP],      gridbox);
res.px_HardCap_woTHP = reshape([R.px_HardCap_woTHP], gridbox);
res.px_HardCap_wTHP  = reshape([R.px_HardCap_wTHP],  gridbox);
res.px_SoftCap_woTHP = reshape([R.px_SoftCap_woTHP], gridbox);
res.px_SoftCap_wTHP  = reshape([R.px_SoftCap_wTHP],  gridbox);
res.px_EVM_woTHP     = reshape([R.px_EVM_woTHP],     gridbox);
res.px_EVM_wTHP      = reshape([R.px_EVM_wTHP],      gridbox);
res.px_EVMref_woTHP  = reshape([R.px_EVMref_woTHP],  gridbox);
res.px_EVMref_wTHP   = reshape([R.px_EVMref_wTHP],   gridbox);

% ★S2P route
res.s2p_BER_woTHP     = reshape([R.s2p_BER_woTHP],     gridbox);
res.s2p_BER_wTHP      = reshape([R.s2p_BER_wTHP],      gridbox);
res.s2p_HardCap_woTHP = reshape([R.s2p_HardCap_woTHP], gridbox);
res.s2p_HardCap_wTHP  = reshape([R.s2p_HardCap_wTHP],  gridbox);
res.s2p_SoftCap_woTHP = reshape([R.s2p_SoftCap_woTHP], gridbox);
res.s2p_SoftCap_wTHP  = reshape([R.s2p_SoftCap_wTHP],  gridbox);
res.s2p_EVM_woTHP     = reshape([R.s2p_EVM_woTHP],     gridbox);
res.s2p_EVM_wTHP      = reshape([R.s2p_EVM_wTHP],      gridbox);
res.s2p_EVMref_woTHP  = reshape([R.s2p_EVMref_woTHP],  gridbox);
res.s2p_EVMref_wTHP   = reshape([R.s2p_EVMref_wTHP],   gridbox);

paramSweep = baseParam;
paramSweep.cutoff_coeff = cutoff_list;
paramSweep.SNRdB        = snr_list;
paramSweep.first_F      = firstF_list;
paramSweep.bmax_initial = bmax_list;
paramSweep.MODNUM       = mod_list;
paramSweep.gridbox      = gridbox;
paramSweep.meshN        = meshN;
paramSweep.winBase      = winBase;
paramSweep.usePartialExp = usePartialExp;

% ★S2P メタ情報も保存
paramSweep.useS2P   = useS2P;
paramSweep.s2pFile  = s2pFile;
paramSweep.Zs = Zs; paramSweep.Zl = Zl;
paramSweep.optS2P = optS2P;
paramSweep.infoS2P_by_cutoff = info_s2p_by_cutoff;  % Bdev 等（軽い情報だけの想定）

outfile = "Sweep3routes_" + num2str(baseParam.TXD_N) + "sym.mat";
save(outfile, 'paramSweep', 'res', '-v7.3');
fprintf('Saved: %s\n', outfile);
