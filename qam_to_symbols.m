function [sym_idx, const] = qam_to_symbols(signal, M)
% qam_to_symbols
%   QAM 波形列を最近傍コンステレーション点に割り当てる。
% 入力:
%   signal : 複素シンボル列
%   M      : QAM次数（正方形）
% 出力:
%   sym_idx: 1..M のインデックス
%   const  : 使用したコンステレーション（列ベクトル）

    const = qam_constellation(M);
    sig = signal(:).';
    sym_idx = zeros(size(sig));
    for n = 1:numel(sig)
        [~, k] = min(abs(sig(n) - const));
        sym_idx(n) = k;
    end
end
