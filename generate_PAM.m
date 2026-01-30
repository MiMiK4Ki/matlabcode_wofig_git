function PAM = generate_PAM(M, TXD_N)
    % 変数の初期化
    A = 1; % 振幅のスケーリングファクタ

    % Mが偶数であることを確認
    if mod(M, 2) ~= 0
        error('M must be even');
    end

    % 0からM-1までのランダムな整数を生成
    TXD = randi([0, M-1], [1, TXD_N]);

    % M-PAM変調を実行
    levels = -(M-1):2:(M-1); % M-PAMのレベル
    PAM = zeros(1, numel(TXD));

    for i = 1:M
        PAM(TXD == (i-1)) = levels(i);
    end

    % 振幅のスケーリング
    PAM = A * PAM;
end
