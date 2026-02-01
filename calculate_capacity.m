function [Capacity, a_hat_] = calculate_capacity(foldername, a_hat, tx_sym, SNRdB_value, first_F_value, first_C, suffix, Mdummy)
% calculate_capacity
%   離散入力(tx_sym)・連続出力(a_hat)の相互情報量をヒストグラムで推定
%   - 実数出力: 1Dヒスト
%   - 複素出力: 2Dヒスト（I/Q平面）
%
% 注意:
%   - AWGN仮定(論文式11)ではなく、経験分布に基づく(式9,10の数値近似)。
%   - suffix, SNRdB_value はファイル名互換のため残す。

    %#ok<NASGU> Mdummy  % 旧引数互換（使わない）

    shift = abs(first_C - first_F_value);

    tx_tr = tx_sym(1:end-shift);
    rx_tr = a_hat(shift+1:end);

    Symbolunique = unique(tx_tr);
    q = numel(Symbolunique);

    if q < 2
        Capacity = NaN;
        a_hat_ = [];
        return;
    end

    % -------------------------
    % 実数出力: 1D
    % -------------------------
    if isreal(rx_tr) && isreal(tx_tr)
        BinEdges = linspace(min(rx_tr), max(rx_tr), 150);
        BinWidth = BinEdges(2) - BinEdges(1);

        a_hat_ = struct('symbolunique', cell(1,q), 'histbin', [], 'Values', [], 'BinWidth', BinWidth);

        P = zeros(q, numel(BinEdges)-1);  % p(y|xk) をビン上で保持
        for k = 1:q
            data = rx_tr(tx_tr == Symbolunique(k));
            [counts, edges] = histcounts(data, BinEdges, 'Normalization','pdf');
            centers = (edges(1:end-1) + edges(2:end))/2;

            a_hat_(k).symbolunique = Symbolunique(k);
            a_hat_(k).histbin      = centers;
            a_hat_(k).Values       = counts;

            P(k,:) = counts;
        end

        isum = sum(P, 1);   % Σ_k p(y|xk)

        ksum = 0;
        for k = 1:q
            pk = P(k,:);
            mask = (pk > 0) & (isum > 0);
            ksum = ksum + sum( pk(mask) .* log2(isum(mask)./pk(mask)) ) * BinWidth;
        end

        Capacity = log2(q) - (1/q)*ksum;

    % -------------------------
    % 複素出力: 2D (I/Q)
    % -------------------------
    else
        Nb = 80;  % ビン数（データが多いなら増やせる）
        xr = real(rx_tr); xi = imag(rx_tr);

        % 範囲（外れ値が強いなら quantile で切るのも有効）
        re_min = min(xr); re_max = max(xr);
        im_min = min(xi); im_max = max(xi);

        BinEdgesRe = linspace(re_min, re_max, Nb+1);
        BinEdgesIm = linspace(im_min, im_max, Nb+1);
        dRe = BinEdgesRe(2) - BinEdgesRe(1);
        dIm = BinEdgesIm(2) - BinEdgesIm(1);
        BinArea = dRe * dIm;

        a_hat_ = struct('symbolunique', cell(1,q), 'H', [], 'BinEdgesRe', BinEdgesRe, 'BinEdgesIm', BinEdgesIm, 'BinArea', BinArea);

        P = zeros(Nb, Nb, q);  % p(y|xk) in 2D
        for k = 1:q
            data = rx_tr(tx_tr == Symbolunique(k));
            H = histcounts2(real(data), imag(data), BinEdgesRe, BinEdgesIm, 'Normalization','pdf');

            P(:,:,k) = H;
            a_hat_(k).symbolunique = Symbolunique(k);
            a_hat_(k).H = H;
        end

        isum = sum(P, 3);  % Σ_k p(y|xk)

        ksum = 0;
        for k = 1:q
            pk = P(:,:,k);
            mask = (pk > 0) & (isum > 0);
            ksum = ksum + sum( pk(mask) .* log2(isum(mask)./pk(mask)) ) * BinArea;
        end

        Capacity = log2(q) - (1/q)*ksum;
    end

    % 保存（旧互換）
    filename = fullfile(foldername, sprintf('capacity%s_fF%d_snr%d.mat', suffix, first_F_value, SNRdB_value));
    save(filename, 'Capacity');
end
