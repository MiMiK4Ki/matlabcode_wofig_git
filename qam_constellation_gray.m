function const = qam_constellation_gray(M)
% qam_constellation_gray
%   標準の正方格子QAM（未正規化）を2軸Gray順で並べた const(1..M) を返す。

    L = sqrt(double(M));
    if abs(L - round(L)) > 1e-12
        error('qam_constellation_gray: M must be a perfect square.');
    end
    L = round(L);

    b = log2(double(L));
    if abs(b - round(b)) > 1e-12
        error('qam_constellation_gray: Gray mapping needs sqrt(M) to be power of 2 (M=4,16,64,256,...).');
    end
    b = round(b);

    levels = double(-(L-1):2:(L-1));  % 例: L=8 -> [-7 -5 ... 7]

    const = zeros(M, 1);
    for q = uint32(0):(uint32(L)-1)
        gQ = bitxor(q, bitshift(q, -1));  % gray
        for i = uint32(0):(uint32(L)-1)
            gI = bitxor(i, bitshift(i, -1));

            g = bitshift(gQ, b) + gI;
            const(double(g)+1) = levels(double(i)+1) + 1j*levels(double(q)+1);
        end
    end
end
