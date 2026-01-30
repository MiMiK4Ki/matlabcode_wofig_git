function D = io_discretize(IO, coeffs, ref_sym)
% io_discretize
%   連続波形 IO.y_wf から、評価/LS推定に使う「受信シンボル中心サンプル列」を作る
%   旧規約に合わせて idx = c_index + (first_C + k)*NoSpS を採用
%
% 入力:
%   IO.y_wf, IO.NoSpS, IO.c_index
%   coeffs.first_C
%   ref_sym : 参照送信シンボル列（長さK決定用）
%
% 出力:
%   D.y_sym   : 受信シンボル中心サンプル列（raw, まだNormRefで割っていない）
%   D.idx     : y_wf 内のサンプル位置
%   D.K       : 使用シンボル数
%   D.c_index : 使用した c_index

y      = IO.y_wf(:).';
NoSpS  = IO.NoSpS;
c_index = IO.c_index;
first_C = coeffs.first_C;

% idx の最後が numel(y) を超えない範囲で K を決める（first_C を考慮）
% max idx = c_index + (first_C + (K-1))*NoSpS <= numel(y)
Kmax = floor((numel(y) - c_index)/NoSpS - first_C) + 1;

Kref = numel(ref_sym);
K = min(Kref, Kmax);

idx = c_index + (first_C + (0:K-1))*NoSpS;

D = struct();
D.y_sym   = y(idx);
D.idx     = idx;
D.K       = K;
D.c_index = c_index;
end
