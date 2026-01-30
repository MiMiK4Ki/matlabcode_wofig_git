function [bh_square, cross_blbk] = calc_bh_crossterm(bk, h_array, deltat,N)
% calc_bh_crossterm
%   bk : シンボル列 (1×N)
%   h_array : 長さNのセル配列 h_array{k} = 波形データ(ベクトル)
%            例: h_array{k}(m) ~ h(t_m - kT_c)
%   deltat : サンプリング刻み [sec]
%
% 出力:
%   bh_square(k) = ∫ (bk·h_k(t))^2 dt (数値和)
%   cross_blbk(k,l) = ∫(bk·bl·h_k(t)·h_l(t)) dt,  k≠l


    bh_square   = zeros(1,N);
    cross_blbk = zeros(N,N);

    for kSym = 1:N
        % (bk·h_k)^2 を積分
        yk = bk(kSym)* h_array{kSym};         % 波形(ベクトル)
        bh_square(kSym) = sum( (yk).^2 ) * deltat;
    end

    % k≠l のクロスターム
    for kSym = 1:N-1
        for lSym = kSym+1:N
            yk = bk(kSym)* h_array{kSym};
            yl = bk(lSym)* h_array{lSym};
            cross_val = sum( yk .* yl ) * deltat; 
            cross_blbk(kSym,lSym) = cross_val;
            cross_blbk(lSym,kSym) = cross_val;  %対称に
        end
    end
end
