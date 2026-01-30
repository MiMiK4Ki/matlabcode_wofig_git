function saveFoldedGaussian(foldername, a_hat_, fitResults, first_F_value, SNRdB_value, suffix)
% SAVEFOLDEDGAUSSIAN Save folded Gaussian data without plotting figures.

    for i = 1:numel(a_hat_)
        sym = a_hat_(i).symbolunique;
        data = a_hat_(i).data;

        centers = a_hat_(i).histbin;
        bw = a_hat_(i).BinWidth;
        edges = [centers - bw/2, centers(end) + bw/2];

        idx = find([fitResults.symbol] == sym, 1);
        if isempty(idx)
            warning('Symbol %d のフィット結果が見つかりません', sym);
            continue;
        end
        xFold = fitResults(idx).xFold;
        yFold = fitResults(idx).yFold;

        filename = sprintf('foldedFit%s_sym%d_fF%d_snr%d.mat', suffix, sym, first_F_value, SNRdB_value);
        save(fullfile(foldername, filename), 'data', 'edges', 'xFold', 'yFold');

    end
end
