function bits = qam_to_bits(sym_idx, M)
% qam_to_bits
%   QAMシンボルインデックス (1..M) をビット列へ変換。
%   現状は自然二進数 (left-msb) の対応。

    if mod(sqrt(M), 1) ~= 0
        error('M must be a perfect square for square QAM.');
    end

    idx0 = sym_idx(:) - 1;
    bits_per_symbol = log2(M);
    bits = de2bi_local(idx0, bits_per_symbol, 'left-msb');
    bits = reshape(bits.', 1, []);
end

function bits = de2bi_local(d, n, varargin)
    bits = zeros(length(d), n);
    for i = n:-1:1
        bits(:,i) = mod(d, 2);
        d = floor(d/2);
    end
    if nargin == 3 && strcmp(varargin{1}, 'left-msb')
        bits = fliplr(bits);
    end
end
