function bits = qam_to_bits(sym_idx, M)
% qam_to_bits (Gray mapping)
%   sym_idx は qam_to_symbols の出力（Gray順インデックス）を想定。
%   idx0 = sym_idx-1 をそのままビット化すれば Grayビット列になる。

    bits_per_symbol = log2(M);
    idx0 = uint32(sym_idx(:) - 1);

    bits_mat = zeros(numel(idx0), bits_per_symbol, 'uint8');
    for b = 1:bits_per_symbol
        shift = bits_per_symbol - b;
        bits_mat(:, b) = uint8(bitand(bitshift(idx0, -shift), 1));
    end

    bits = reshape(bits_mat.', 1, []);
end
