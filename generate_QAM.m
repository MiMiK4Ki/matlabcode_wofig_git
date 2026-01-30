function QAM = generate_QAM(M, TXD_N)
% generate_QAM
%   正方形 QAM (M=4,16,64...) のシンボル列を生成。
%   既定では平均電力の正規化は行わず、格子点そのままを返す。

    if mod(sqrt(M), 1) ~= 0
        error('M must be a perfect square for square QAM.');
    end

    const = qam_constellation(M);
    idx = randi([1, M], [1, TXD_N]);
    QAM = const(idx);
end
