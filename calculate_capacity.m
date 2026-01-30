function [Capacity , a_hat_] = calculate_capacity(foldername,a_hat,PAM_signal,SNRdB_value,first_F_value,first_C,suffix,M)
% CALCULATE_CAPACITY Evaluate capacity without plotting graphs.

PAM_truncated = PAM_signal(1:end-abs(first_C-first_F_value));
a_hat_truncated = a_hat(abs(first_C-first_F_value)+1:end);

Symbolunique = unique(PAM_signal); % PAM symbol values
BinEdges = linspace(min(a_hat_truncated), max(a_hat_truncated), 150);
BinWidth = BinEdges(2) - BinEdges(1);

a_hat_ = struct('symbolunique', cell(1, numel(Symbolunique)), ...
                'data', [], 'histbin', [], 'Values', [], 'BinWidth', []);

for uniqueindex = 1:numel(Symbolunique)
    a_hat_(uniqueindex).symbolunique = Symbolunique(uniqueindex);
    data = a_hat_truncated(PAM_truncated == Symbolunique(uniqueindex));
    a_hat_(uniqueindex).data = data;

    [counts, edges] = histcounts(data, BinEdges, 'Normalization','pdf');
    centers = (edges(1:end-1) + edges(2:end))/2;

    a_hat_(uniqueindex).histbin = centers;
    a_hat_(uniqueindex).Values  = counts;
    a_hat_(uniqueindex).BinWidth = edges(2) - edges(1);
end

isum = zeros(1,numel(a_hat_(1).Values));
for idx = 1:numel(Symbolunique)
    isum = isum + a_hat_(idx).Values;
end

ksum = 0;
for idx = 1:numel(Symbolunique)
    ksum = ksum + nansum( a_hat_(idx).Values .* log2(isum./a_hat_(idx).Values) ) .* BinWidth;
end

Capacity = log2(numel(Symbolunique)) - 1/numel(Symbolunique) * ksum;

% 保存
filename = fullfile(foldername, sprintf('capacity%s_fF%d_snr%d.mat', suffix, first_F_value, SNRdB_value));
save(filename, 'Capacity');

end