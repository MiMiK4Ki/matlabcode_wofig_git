function TXD = pam_to_bits(PAM, Modnum, A)
    if nargin < 3
        A = 1; % スケーリングファクタが指定されていない場合、デフォルト値1を使用
    end

    % Mが偶数であることを確認
    if mod(Modnum, 2) ~= 0
        error('M must be even');
    end

    % M-PAMのレベル
    levels = -(Modnum-1):2:(Modnum-1);
    
    % スケーリングを元に戻す
    PAM = PAM / A;
    
    % 振幅値からシンボルインデックスを特定
    [~, indices] = min(abs(PAM - levels'), [], 1);
    indices = indices' - 1; % 0からM-1までのインデックスに変換
    
    % シンボルインデックスからビット列を生成
    bits_per_symbol = log2(Modnum);
    TXD = de2bi(indices, bits_per_symbol, 'left-msb');
    TXD = reshape(TXD',1,[]); % 1次元配列に変換
end

function bits = de2bi(d, n, varargin)
    % デシマル値をビット列に変換する関数
    bits = zeros(length(d), n);
    for i = n:-1:1
        bits(:,i) = mod(d, 2);
        d = floor(d/2);
    end
    if nargin == 3 && strcmp(varargin{1}, 'left-msb')
        bits = fliplr(bits);
    end
end
