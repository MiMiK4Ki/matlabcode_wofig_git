function [sym_idx, levels] = pam_to_symbols(signal, M, normalize_factor)
% pam_to_symbols  M-PAM の振幅系列をシンボル番号 (1:M) に戻す
%  入力:
%    signal           - 入力波形列 (1×N または N×1)
%    M                - PAM のモダリング数 (例: 16)
%    normalize_factor - (オプション) 信号正規化係数。省略時は 1.
%  出力:
%    sym_idx - 受信シンボルのインデックス列 (1～M)
%    levels  - 各サンプルの再マッピングされた PAM レベル

    if nargin < 3
        normalize_factor = 1;
    end

    % M-PAM の理想振幅レベル: 例えば M=4 なら [-3,-1,+1,+3]
    levels_ideal = (-(M-1):2:(M-1)) * normalize_factor;  

    % 各サンプルごとに「最も近い」振幅レベルを探し、シンボル番号を決定
    N = numel(signal);
    sym_idx = zeros(size(signal));
    levels  = zeros(size(signal));
    for n = 1:N
        [~, k] = min(abs(signal(n) - levels_ideal));
        sym_idx(n) = k;
        levels(n)  = levels_ideal(k);
    end
end
