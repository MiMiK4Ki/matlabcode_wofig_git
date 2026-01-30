function const = qam_constellation(M)
% qam_constellation
%   正方形 QAM のコンステレーション点を返す（列ベクトル）。
%   並び順は I 軸の小→大、Q 軸の小→大（meshgridの行優先）。

    if mod(sqrt(M), 1) ~= 0
        error('M must be a perfect square for square QAM.');
    end

    sqrtM = sqrt(M);
    levels = -(sqrtM-1):2:(sqrtM-1);
    [I, Q] = meshgrid(levels, levels);
    const = I(:) + 1j * Q(:);
end
