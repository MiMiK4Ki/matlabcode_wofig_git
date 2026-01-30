%% compute_transition_and_hardcapacity.m
% このファイルには遷移確率行列 r_ik と Hard Capacity を計算する関数を定義します。

function [r, C_hard] = compute_transition_and_hardcapacity(TX_sym, RX_sym, M)
% compute_transition_and_hardcapacity  遷移確率行列 r(i,k) と Hard Capacity を計算
%   入力:
%     TX_sym : 送信シンボル列 (1×N または N×1, 値は 1 から M)
%     RX_sym : 受信シンボル列 (同上)
%     M      : シンボル数
%   出力:
%     r      : M×M 遷移確率行列, r(i,k) = P(rxn=k | tx=i)
%     C_hard : Hard Capacity

    % 1) カウント行列の構築
    counts = zeros(M, M);
    for i = 1:M
        idx = find(TX_sym == i);
        if isempty(idx)
            continue;
        end
        for k = 1:M
            counts(i, k) = sum(RX_sym(idx) == k);
        end
    end

    % 2) 遷移確率行列 r(i,k) の正規化 (行ごとに合計1)
    row_sums = sum(counts, 2);
    r = zeros(M, M);
    for i = 1:M
        if row_sums(i) > 0
            r(i, :) = counts(i, :) / row_sums(i);
        end
    end

    % 3) Hard Capacity の計算
    % C = log2(M) + (1/M) * sum_{i=1}^M sum_{k=1}^M r(i,k) * log2( r(i,k) / sum_j r(j,k) )
    col_sums = sum(r, 1);
    term = 0;
    for i = 1:M
        for k = 1:M
            if r(i,k) > 0 && col_sums(k) > 0
                term = term + r(i,k) * log2( r(i,k) / col_sums(k) );
            end
        end
    end
    C_hard = log2(M) + term / M;
end

%% main.m に組み込む例
% （BER計算部分の直後に追加）
% [TX_sym, ~] = pam_to_symbols(PAM_signal, Modnum_value);
% [RX_sym, ~] = pam_to_symbols(a_hat_wTHP, Modnum_value);
% [r_ik, Hard_capacity_temp] = compute_transition_and_hardcapacity(TX_sym, RX_sym, Modnum_value);
% iterationRes.r_ik = r_ik;
% iterationRes.Hard_capacity = Hard_capacity_temp;
