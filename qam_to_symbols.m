function [sym_idx, const] = qam_to_symbols(signal, M, ref_sym)
% qam_to_symbols (Gray mapping)
%   ref_sym を渡した場合、ref_sym から推定したレベルで
%   Gray順に const を構築して最近傍割当する。
%
% 入力:
%   signal : 複素シンボル列
%   M      : QAM次数（正方形）
%   ref_sym: (任意) 参照送信シンボル列（スケール合わせ用）
%
% 出力:
%   sym_idx: 1..M のインデックス（この idx-1 が Grayラベル）
%   const  : Gray順コンステレーション（列ベクトル）

    if mod(sqrt(M),1) ~= 0
        error('qam_to_symbols: M must be a perfect square (square QAM).');
    end

    if nargin >= 3 && ~isempty(ref_sym)
        const = qam_constellation_gray_from_ref(ref_sym, M);
    else
        % ref無しの場合は標準レベル（未正規化）で作る
        const = qam_constellation_gray(M);
    end

    sig = signal(:);
    % 最近傍（N×M行列を作るが M<=256 くらいなら現実的）
    [~, sym_idx] = min(abs(sig - const.'), [], 2);
    sym_idx = sym_idx(:).';  % row
end
