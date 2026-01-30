function coeffs = estimate_channel_ls(IO_or_rx, x_sym, win)
% estimate_channel_ls
%   連続受信波形(IO.y_wf) + 送信シンボル列(x_sym) から
%   symbol-spaced のチャネル係数を LS で推定して coeffs を返す。
%
% 重要:
%   - 離散化(サンプリング)は evaluate_from_io と同じ規約に揃えるため io_discretize を使う
%   - ここで推定するのは「多めの causal taps」をまず推定し、
%     その最大振幅タップを m=0 とみなして first_C:C_max / first_F:F_max を切り出す
%
% 入力:
%   IO_or_rx :
%     - struct( y_wf, NoSpS, c_index ) 連続波形
%     - または rx_sym (受信シンボル列)
%   x_sym : 送信シンボル列（row/colどちらでもOK）
%   win : struct('first_C',-2,'C_max',15,'first_F',0,'F_max',15) など
%
% 出力 coeffs (chan_cont_to_discrete と同じ形):
%   coeffs.ChannelCoeff  : [first_C:C_max] / NormRef
%   coeffs.FilterCoeff   : [first_F:F_max] / NormRef, 先頭は -1 して旧IIR規約に合わせる
%   coeffs.first_C, coeffs.first_F, coeffs.NormRef
%   （デバッグ用に coeffs.h_est_full, coeffs.k_peak, coeffs.m_axis_full も付ける）

% ---- (1) 連続→離散（evaluate_from_io と同じ規約）----
x_sym = x_sym(:).';  % row
if isstruct(IO_or_rx) && isfield(IO_or_rx,'y_wf')
    tmp = struct('first_C', win.first_C);     % io_discretize が必要とする最小フィールド
    D   = io_discretize(IO_or_rx, tmp, x_sym);
    y_sym = D.y_sym(:).';                     % row (未正規化)
    Nsym = D.K;
else
    y_sym = IO_or_rx(:).';
    Nsym = min(numel(x_sym), numel(y_sym));
end

% 長さ合わせ（安全側）
x = x_sym(1:Nsym);
y = y_sym(1:Nsym);

% ---- (2) LS 用の「多め taps」長を決める（かなり素朴）----
% ここは「多めに取る」ための固定 extra
extra = 20;
maxTap = max(win.C_max, win.F_max);
Ntaps = abs(win.first_C) + maxTap + 1 + extra;

% ---- (3) 畳み込み行列 X を作って LS（causal taps 0..Ntaps-1）----
% y(n) ≈ sum_{k=0}^{Ntaps-1} h(k) * x(n-k)
X = zeros(Nsym, Ntaps);
for n = 1:Nsym
    kk = min(Ntaps, n);                 % n-kk+1 >= 1 の範囲だけ埋める
    X(n, 1:kk) = x(n:-1:n-kk+1);
end

% 過渡影響を避けるために中央だけ使う（複雑な条件分岐はしない）
skip = Ntaps;                           % "少し大きめ"のマージン
rows = (skip+1):(Nsym-skip);

h_est = X(rows,:) \ (y(rows).');        % column (Ntaps×1)

% ---- (4) 最大振幅タップを m=0 とみなして、そこから切り抜く ----
[~, kpk] = max(abs(h_est));             % peak index in 1..Ntaps

mG = win.first_C:win.C_max;
hG = zeros(1, numel(mG));
for ii = 1:numel(mG)
    k = kpk + mG(ii);                   % m=0 → kpk
    if (k >= 1) && (k <= Ntaps)
        hG(ii) = h_est(k);
    end
end

mF = win.first_F:win.F_max;
hF = zeros(1, numel(mF));
for ii = 1:numel(mF)
    k = kpk + mF(ii);
    if (k >= 1) && (k <= Ntaps)
        hF(ii) = h_est(k);
    end
end

% ---- (5) 旧規約どおり正規化して coeffs を作る ----
NormRef = hF(1);                         % m = first_F のタップ
ChannelCoeff = hG / NormRef;
FilterCoeff  = hF / NormRef;
FilterCoeff(1) = FilterCoeff(1) - 1;     % 旧IIR/THP の規約

coeffs = struct();
coeffs.ChannelCoeff = ChannelCoeff;
coeffs.FilterCoeff  = FilterCoeff;
coeffs.first_C      = win.first_C;
coeffs.first_F      = win.first_F;
coeffs.NormRef      = NormRef;

% ---- デバッグ用（必要なければ消してOK）----
coeffs.h_est_full  = h_est;              % 0..Ntaps-1 の causal taps（未正規化）
coeffs.k_peak      = kpk;                % peak の index
coeffs.m_axis_full = (0:Ntaps-1) - (kpk-1); % peak を 0 とした相対m軸
end
