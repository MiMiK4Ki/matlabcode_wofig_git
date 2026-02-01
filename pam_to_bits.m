function bits = pam_to_bits(PAM, Modnum, A, mapping)
% pam_to_bits
%   PAMの振幅値列 -> ビット列（1×(N*log2(Modnum))）
%
% 入力:
%   PAM     : PAMレベル列（実数想定。複素でも abs で最近傍に落ちるが、基本は実数）
%   Modnum  : PAM次数（通常 2^k）
%   A       : スケーリング係数（省略時 1）
%   mapping : "bin"（自然2進） or "gray"（Gray符号） 省略時 "bin"
%
% 出力:
%   bits : left-msb で並んだビット列（row）

    if nargin < 3 || isempty(A)
        A = 1;
    end
    if nargin < 4 || isempty(mapping)
        mapping = "bin";
    end
    mapping = lower(string(mapping));

    if mod(Modnum, 2) ~= 0
        error('pam_to_bits: Modnum must be even.');
    end

    bps = log2(Modnum);
    if abs(bps - round(bps)) > 1e-12
        error('pam_to_bits: Modnum must be a power of 2 to map to bits.');
    end
    bps = round(bps);

    % M-PAM levels（未正規化）
    levels = double(-(Modnum-1):2:(Modnum-1));   % 例: 8PAM -> [-7 -5 ... 7]

    % スケールを戻す
    x = PAM(:) / A;     % column

    % 最近傍レベル -> 0..M-1 のインデックス
    [~, idx] = min(abs(x - levels(:).'), [], 2);   % idx: 1..M
    idx0 = uint32(idx - 1);                        % 0..M-1

    % Gray符号化（任意）
    if mapping == "gray"
        idx0 = bitxor(idx0, bitshift(idx0, -1));
    elseif mapping ~= "bin"
        error('pam_to_bits: mapping must be "bin" or "gray".');
    end

    % left-msb のビット列へ（bitgetはLSB=1）
    bits_mat = zeros(numel(idx0), bps, 'uint8');
    for b = 1:bps
        pos = bps - b + 1;            % MSBから取る
        bits_mat(:, b) = uint8(bitget(idx0, pos));
    end

    bits = reshape(bits_mat.', 1, []);
end
