function ChanOut = chan_cont_to_discrete(param, ChanIn, win)
% 連続IR h_rx から、旧規約( first_C:C_max / first_F:F_max )で離散タップを抽出
% 入力:
%   param : struct（少なくとも NoSpS, first_F を持つのが望ましい）
%   ChanIn: struct（h_rx: *dt 済み 合成IR, dt）
%   win   : struct('first_C',-2,'C_max',15,'first_F',param.first_F,'F_max',15)
% 出力:
%   ChanOut = ChanIn に以下を追加した構造体
%     .ChannelCoeff  : [first_C:C_max] の離散タップ / NormRef
%     .FilterCoeff   : [first_F:F_max] の離散タップ / NormRef, 先頭から 1 を引く
%     .first_C, .first_F, .NormRef

    arguments
        param struct
        ChanIn struct
        win struct = struct('first_C',-2,'C_max',15,'first_F',param.first_F,'F_max',15)
    end

    h_rx = ChanIn.h_rx(:).';
    NoSpS = param.NoSpS;

    % 主タップ（シンボル中心）
    [~, c_index] = max(abs(h_rx));

    % ラグ → サンプル位置
    mG   = win.first_C:win.C_max;
    mF   = win.first_F:win.F_max;
    idxG = c_index + mG*NoSpS;
    idxF = c_index + mF*NoSpS;

    % 範囲外は 0 扱い
    hG = zeros(size(idxG)); okG = idxG>=1 & idxG<=numel(h_rx); hG(okG) = h_rx(idxG(okG));
    hF = zeros(size(idxF)); okF = idxF>=1 & idxF<=numel(h_rx); hF(okF) = h_rx(idxF(okF));

    NormRef      = hF(1);                 % 旧来の正規化基準（m = first_F）
    ChannelCoeff = hG / NormRef;          % [first_C:C_max]
    FilterCoeff  = hF / NormRef;          % [first_F:F_max]
    FilterCoeff(1) = FilterCoeff(1) - 1;  % ZF-THP/IIR 形に合わせる（旧来どおり）

    ChanOut              = ChanIn;
    ChanOut.ChannelCoeff = ChannelCoeff;
    ChanOut.FilterCoeff  = FilterCoeff;
    ChanOut.first_C      = win.first_C;
    ChanOut.first_F      = win.first_F;
    ChanOut.NormRef      = NormRef;
end
