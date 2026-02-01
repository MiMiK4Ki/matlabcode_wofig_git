function const = qam_constellation_gray_from_ref(ref_sym, M)
% qam_constellation_gray_from_ref
%   ref_sym の星座点（スケール込み）から I/Q レベルを推定し、
%   2軸Gray順に const(1..M) を構築する。

    L = sqrt(double(M));
    if abs(L - round(L)) > 1e-12
        error('qam_constellation_gray_from_ref: M must be a perfect square.');
    end
    L = round(L);

    b = log2(double(L));
    if abs(b - round(b)) > 1e-12
        error('qam_constellation_gray_from_ref: Gray mapping needs sqrt(M) to be power of 2 (M=4,16,64,256,...).');
    end
    b = round(b);

    pts = unique(ref_sym(:));
    lvI = sort(unique(real(pts)));
    lvQ = sort(unique(imag(pts)));

    % ref_sym が全レベルを含まない場合のフォールバック（安全策）
    if numel(lvI) ~= L || numel(lvQ) ~= L
        const = qam_constellation_gray(M);  % 標準格子（未正規化）へ
        return;
    end

    const = zeros(M, 1);

    for q = uint32(0):(uint32(L)-1)
        gQ = bitxor(q, bitshift(q, -1));  % gray(q)
        for i = uint32(0):(uint32(L)-1)
            gI = bitxor(i, bitshift(i, -1));  % gray(i)

            g = bitshift(gQ, b) + gI;  % 0..M-1
            const(double(g)+1) = lvI(double(i)+1) + 1j*lvQ(double(q)+1);
        end
    end
end
