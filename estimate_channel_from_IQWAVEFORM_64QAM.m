%% estimate_channel_from_IQWAVEFORM_64QAM.m
% 入力シンボル(tx) → 出力シンボル(rx) のFIRチャネル h[m] をLS推定
% 出力シンボルは「(Tx側)sRRC→実験施設→(Rx側)sRRC」後のシンボル列を想定

clear; close all; clc;

%% ===== ユーザ設定 =====
matFile = 'IQWAVEFORM_64QAM.mat';   % ←パスを適宜変更

mRange  = -3:15;     % 推定するタップのラグ（例：-3:15）
skipSym = 200;       % 推定に使わない前後シンボル数（フィルタ過渡回避）
doPhaseCorrection    = true;  % 共通位相回転補正するか
phaseThDeg           = 1.0;   % |推定位相|がこの値[deg]を超えたら補正

%% ===== ロード =====
S = load(matFile);
fprintf('[info] Loaded: %s\n', matFile);

%% ===== 変数の自動選択 =====
% 入力シンボル（優先順）
tx = pickVar(S, {'tx','tx_upd','inp_symb','Inp_symb','txSymbols','TxSymbols'}, true);

if isempty(tx)
    error('TX symbols not found. "tx" (or tx_upd / inp_symb) が .mat に見つかりませんでした。');
end
tx = tx(:);

% 受信後シンボル列が既にあるならそれを優先
rxSym = pickVar(S, {'rxSymbols_5','rxSymbols','rx_symb_dadj_p1','rx_symb','rxsym'}, true);

% 無ければ、受信波形を拾って sRRC(受信) + ダウンサンプルでシンボル化
if isempty(rxSym)
    rxWave = pickVar(S, {'rx_norm_out1','rx_norm1','xout1','iqdata_vsa1','xout','xout_filt'}, true);
    if isempty(rxWave)
        error(['RX waveform not found. "xout1" (or iqdata_vsa1 など) が .mat に見つかりませんでした。' ...
               newline 'rxSymbols_* が保存されていない場合は受信波形が必要です。']);
    end

    % オーバーサンプリング（nSamps）を取得（無ければ 8）
    if isfield(S,'nSamps')
        Nover = double(S.nSamps);
    elseif isfield(S,'NoverSampling')
        Nover = double(S.NoverSampling);
    else
        Nover = 8;
        fprintf('[warn] nSamps not found -> assume NoverSampling = %d\n', Nover);
    end

    rxSym = AnalyzeQAMWaveform(rxWave(:), Nover);
end
rxSym = rxSym(:);

fprintf('[info] tx symbols = %d, rx symbols = %d\n', numel(tx), numel(rxSym));

%% ===== シンボルレートで遅延合わせ（finddelay）=====
d = local_finddelay(tx, rxSym);
fprintf('[info] Symbol-rate delay (rx relative to tx) = %d [symbols]\n', d);

[txA, rxA] = align_by_delay(tx, rxSym, d);
K = min(numel(txA), numel(rxA));
txA = txA(1:K);
rxA = rxA(1:K);

%% ===== 共通位相回転（CPE）補正（必要なら）=====
phiDeg = angle(sum(rxA .* conj(txA))) * 180/pi;
fprintf('[info] Estimated common phase offset = %+0.3f deg\n', phiDeg);

rxP = rxA;
if doPhaseCorrection && abs(phiDeg) > phaseThDeg
    rxP = rxA .* exp(-1j * deg2rad(phiDeg));
    fprintf('[info] Applied phase de-rotation: exp(-j*%0.3f deg)\n', phiDeg);
else
    fprintf('[info] Phase correction skipped (already small or disabled).\n');
end

%% ===== FIRチャネル推定（LS）=====
h = estimate_fir_ls_sym(txA, rxP, mRange, skipSym);   % row vector
fprintf('\n==== Estimated FIR channel h[m] (mRange = %d:%d) ====\n', min(mRange), max(mRange));
fprintf('  m      |h[m]|        angle(deg)\n');
for ii = 1:numel(mRange)
    fprintf('%+4d   % .6e    % .3f\n', mRange(ii), abs(h(ii)), angle(h(ii))*180/pi);
end

%% ===== 推定精度（NMSE）を計算（estimate関数と同じ有効区間）=====
[nmse_dB, rxRef, rxHat] = eval_nmse_on_valid_region(txA, rxP, h, mRange, skipSym);
fprintf('\n[info] NMSE on valid region = %0.2f dB\n', nmse_dB);

%% ===== 可視化 =====
figure('Name','Estimated channel taps |h[m]|');
stem(mRange, real(h), 'filled'); grid on;
xlabel('tap lag m'); ylabel('|h[m]|'); title('|h[m]|');

figure('Name','Estimated channel taps phase');
stem(mRange, angle(h)*180/pi, 'filled'); grid on;
xlabel('tap lag m'); ylabel('phase [deg]'); title('angle(h[m])');

figure('Name','Model check: rxRef vs rxHat (constellation-like)');
plot(real(rxRef), imag(rxRef), '.', 'DisplayName','rx (ref)'); hold on;
plot(real(rxHat), imag(rxHat), '.', 'DisplayName','rx\_hat (Xh)'); grid on; axis equal;
xlabel('I'); ylabel('Q'); legend('Location','best');
title(sprintf('Model check (NMSE = %0.2f dB)', nmse_dB));

%% ===== 周波数応答（任意）=====
Nfft = 2048;
H = fftshift(fft(h, Nfft));
f = linspace(-0.5, 0.5, Nfft); % 正規化周波数（シンボルレート基準）

figure('Name','Estimated channel frequency response');
plot(f, 20*log10(abs(H)+1e-12)); grid on;
xlabel('Normalized frequency (cycles/symbol)'); ylabel('|H| [dB]');
title('Estimated channel frequency response');

%% ===== 必要なら保存 =====
outFile = 'estimated_channel_64QAM.mat';
save(outFile, 'h', 'mRange', 'skipSym', 'phiDeg', 'nmse_dB');
fprintf('[info] Saved: %s\n', outFile);

%% ====== local functions ======
function x = pickVar(S, nameList, mustBeVector)
    x = [];
    for i = 1:numel(nameList)
        nm = nameList{i};
        if isfield(S, nm)
            cand = S.(nm);
            if isempty(cand), continue; end
            if ~isnumeric(cand), continue; end
            if mustBeVector && ~isvector(cand), continue; end
            x = cand;
            fprintf('[info] Using variable "%s"\n', nm);
            return;
        end
    end
end

function d = local_finddelay(x, y)
    % finddelay が無い環境用に簡易フォールバック
    if exist('finddelay','file') == 2
        d = finddelay(x, y);
        return;
    end
    % xcorrで最大点
    [c,lags] = xcorr(y, x); % y relative to x
    [~,idx] = max(abs(c));
    d = lags(idx); % y is delayed by d relative to x（正のときyが遅い扱い）
end

function [xA, yA] = align_by_delay(x, y, d)
    % d >= 0: y が x より d サンプル遅い → y先頭を捨てる
    % d < 0 : x が y より -d 遅い → x先頭を捨てる
    if d >= 0
        yA = y(d+1:end);
        xA = x(1:min(end, numel(yA)));
    else
        xA = x((-d)+1:end);
        yA = y(1:min(end, numel(xA)));
    end
end

function ysym = AnalyzeQAMWaveform(x, NoverSampling)
    % 受信sRRC（RaisedCosineReceiveFilter）+ decimation でシンボル化
    nSym = 32;       % フィルタ長（symbols）
    beta = 0.35;     % roll-off
    FilterGain = 0.3578;

    if exist('comm.RaisedCosineReceiveFilter','class') == 8
        hRxFilter = comm.RaisedCosineReceiveFilter( ...
            'RolloffFactor', beta, ...
            'InputSamplesPerSymbol', NoverSampling, ...
            'FilterSpanInSymbols', nSym, ...
            'DecimationFactor', NoverSampling, ...
            'Gain', 1/FilterGain);

        ysym = step(hRxFilter, x);
    else
        % Comm Toolboxが無い場合のフォールバック（rcosdesign）
        % ここでは sqrt-RRC を使う（近い挙動にする意図）
        hr = rcosdesign(beta, nSym, NoverSampling, 'sqrt') / FilterGain;
        y  = filter(hr, 1, x);
        ysym = y(1:NoverSampling:end);
    end
end

function h = estimate_fir_ls_sym(tx, rx, mRange, skip)
    % rx[k] ≈ Σ_{m in mRange} h[m] * tx[k - m] を最小二乗で推定
    if nargin < 4 || isempty(skip)
        skip = 50;
    end

    tx = tx(:);
    rx = rx(:);

    K = min(numel(tx), numel(rx));
    tx = tx(1:K);
    rx = rx(1:K);

    mRange = mRange(:).';
    m_min = min(mRange);
    m_max = max(mRange);

    k_start = 1 + m_max;
    k_end   = K + m_min;

    k_start = k_start + skip;
    k_end   = k_end   - skip;

    if k_end <= k_start
        error('Not enough samples after skipping. skip=%d is too large for current data length.', skip);
    end

    k = (k_start:k_end).';

    X = zeros(numel(k), numel(mRange));
    for ii = 1:numel(mRange)
        X(:,ii) = tx(k - mRange(ii));
    end

    h = X \ rx(k);
    h = h(:).';
end

function [nmse_dB, rxRef, rxHat] = eval_nmse_on_valid_region(tx, rx, h, mRange, skip)
    tx = tx(:); rx = rx(:);
    K = min(numel(tx), numel(rx));
    tx = tx(1:K); rx = rx(1:K);

    mRange = mRange(:).';
    m_min = min(mRange);
    m_max = max(mRange);

    k_start = 1 + m_max + skip;
    k_end   = K + m_min - skip;

    k = (k_start:k_end).';
    X = zeros(numel(k), numel(mRange));
    for ii = 1:numel(mRange)
        X(:,ii) = tx(k - mRange(ii));
    end

    rxRef = rx(k);
    rxHat = X * h(:);

    e = rxRef - rxHat;
    nmse = (norm(e)^2) / (norm(rxRef)^2 + eps);
    nmse_dB = 10*log10(nmse + eps);
end
