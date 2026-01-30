function [xnorm, ynorm, evm, evmref]  = CalcEVMv2_silent(rx,tx,qamOrder)
% CalcEVMv2_silent
%   実験チーム提供 CalcEVMv2 と同じ計算（RMS正規化 + xcorr遅延補償 + 位相補償）
%   ただしコマンドウィンドウへの表示をしない版。
%
% 入力:
%   rx, tx   : ベクトル（実数/複素どちらでもOK）
%   qamOrder : (任意) evmref 用の QAM 次数。既定=64
%
% 出力:
%   xnorm, ynorm : CalcEVMv2 同様に正規化された列
%   evm          : Tx参照EVM（comm.EVM(ynorm,xnorm)）
%   evmref        : 参照星座（QAM）から推定したEVM（nargout>=4のとき計算）

    if nargin < 3 || isempty(qamOrder)
        qamOrder = 64;
    end

    hEVM = comm.EVM;

    StartUpDelay    = 1;
    SamplesPerSymbol = 1;

    x = rx(StartUpDelay+1:SamplesPerSymbol:end);
    y = tx(StartUpDelay+1:SamplesPerSymbol:end);

    c = xcorr(x,y);
    [XcorrPeak, XcorrMaxInd] = max(c);

    % CalcEVMv2 と同じ式（※ mean(numel(x),numel(y)) は numel(x) と同義になる点は仕様踏襲）
    DelayInd = XcorrMaxInd - (mean(numel(x), numel(y))) + 1;
    if DelayInd < 1
        DelayInd = 1;
    end

    xDelayed   = x(DelayInd:1:end);
    yCancelled = y(1:1:(end-DelayInd+1));

    SelForNInd = find(abs(yCancelled) <= 100000);
    yForRot    = yCancelled(SelForNInd);
    xForRot    = xDelayed(SelForNInd);

    angleMean = angle(XcorrPeak);

    % RMS正規化（v2）
    xnorm = xDelayed./(rms(xForRot)).*exp(-1i*angleMean);
    ynorm = yCancelled./(rms(yForRot));

    evm = hEVM(ynorm,xnorm);

    % 参照星座EVM（要求されたときだけ）
    evmref = NaN;
    if nargout >= 4
        hEVM2 = comm.EVM;
        hEVM2.ReferenceSignalSource = "Estimated from reference constellation";
        hEVM2.Normalization = 'Average reference signal power';

        refconst = qammod(0:(qamOrder-1), qamOrder);
        refconst = refconst./rms(refconst)*rms(abs(xnorm));   % CalcEVMv2 と同じスケーリング
        hEVM2.ReferenceConstellation = refconst;

        evmref = hEVM2(xnorm);
    end
end
