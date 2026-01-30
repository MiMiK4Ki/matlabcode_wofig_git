function data = calculate_dsi(dk, sk, ik, P_value, epsilon, varargin)
    % Normalized_Coeffの設定
    if nargin > 5
        Normalized_Coeff = varargin{1} ./ sqrt(P_value);
    else
        Normalized_Coeff = 1 ./ sqrt(P_value);
    end
    
    % dkの最も近い2つの数値とその距離、および種類を数える
    [min_distance, ~, unique_count] = min_distance_closest_pair(dk, epsilon);
    
    % SignalDの計算
    data.SignalD = min_distance .*abs(Normalized_Coeff);
    data.ModNum = unique_count;

    % 各種分散の計算
    data.dk = var(dk .* Normalized_Coeff);
    data.sk = var(sk .* Normalized_Coeff);
    data.ik = var(ik .* Normalized_Coeff);

    data.absdk = var(abs(dk) .* Normalized_Coeff);
    data.abssk = var(abs(sk) .* Normalized_Coeff);
    data.absik = var(abs(ik) .* Normalized_Coeff);
end


function [min_distance, closest_pair, unique_count] = min_distance_closest_pair(dk, epsilon)
    % dkを昇順にソート
    dk_sorted = sort(unique(dk));  % 重複を排除して昇順にソート

    % 隣接する要素間の距離を計算
    distances = diff(dk_sorted);

    % 極微小な差を無視するインデックスを取得
    valid_indices = find(distances > epsilon);

    if isempty(valid_indices)
        error('有効な距離がありません。すべての要素がほぼ同じです。');
    end

    % 最小の距離を特定
    [min_distance, min_index] = min(distances(valid_indices));

    % 元のインデックスを取得
    original_min_index = valid_indices(min_index);

    % 最も近い2つの数値を取得
    closest_pair = dk_sorted(original_min_index:original_min_index + 1);

    % 近すぎる値を同じ値とみなして種類を数える
    unique_dk = dk_sorted(1);
    for i = 2:length(dk_sorted)
        if dk_sorted(i) - unique_dk(end) > epsilon
            unique_dk = [unique_dk, dk_sorted(i)];
        end
    end
    unique_count = length(unique_dk);
end


