function h = srrc_filter_wosym(t, r, Tc, B)
    % t: 時間ベクトル
    % r: ロールオフファクター
    % Tc: シンボル期間
    % B: 帯域幅調整係数

    h = zeros(size(t));
    eps_t    = 1e-12*Tc;
    % 中央の値（t == 0）の計算
    idx_zero = abs(t) < eps_t;
    h(idx_zero) = (1 - r) + 4 * r / pi;

    % 特殊ケース（t == ±Tc/(4r)）の計算
    idx_special = abs(t - Tc/(4*r)) < eps_t | abs(t + Tc/(4*r)) < eps_t;

    h(idx_special) = r / sqrt(2) * ((1 + 2/pi) * sin(pi/(4*r)) + (1 - 2/pi) * cos(pi/(4*r)));

    % その他のケースの計算
    idx_else = ~(idx_zero | idx_special);
    t_else = t(idx_else);
    h(idx_else) = (sin(pi * t_else / Tc * (1 - r)) + 4 * r * t_else / Tc .* cos(pi * t_else / Tc * (1 + r))) ...
                  ./ (pi * t_else / Tc .* (1 - 16 * r^2 * t_else.^2 / Tc^2));

    % フィルタのスケーリング
    h = h * B / Tc;
end
