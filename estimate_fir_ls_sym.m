function h = estimate_fir_ls_sym(tx, rx, mRange)
% estimate_fir_ls_sym
%   rx[k] ≈ Σ_{m in mRange} h[m] * tx[k - m] を最小二乗で推定
%
% 入力:
%   tx, rx  : 同期済みのシンボル列（同じ長さ推奨）
%   mRange  : 推定するタップのラグ（例: -3:15）
%
% 出力:
%   h : mRange順の推定係数

tx = tx(:);
rx = rx(:);

K = min(numel(tx), numel(rx));
tx = tx(1:K);
rx = rx(1:K);

mRange = mRange(:).';
m_min = min(mRange);
m_max = max(mRange);

% 有効な k 範囲（全部の tx(k-m) が範囲内になるところ）
k_start = 1 + m_max;
k_end   = K + m_min;

% 簡単に安全マージン（トランジェント回避のため）
skip = 50;
k_start = k_start + skip;
k_end   = k_end   - skip;

k = (k_start:k_end).';

X = zeros(numel(k), numel(mRange));
for ii = 1:numel(mRange)
    X(:,ii) = tx(k - mRange(ii));
end

% 最小二乗（MATLABの \ が内部でQR/SVDで解く）
h = X \ rx(k);

h = h(:).';  % rowで返す
end
