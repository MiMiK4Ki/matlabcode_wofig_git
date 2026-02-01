function TX = thp_prepare_tx(param, coeffs, x_sym)
% thp_prepare_tx  実験送信用の THP 列を生成（IIRoutput を外で呼ぶ）
% 入力:
%   param  : struct（MODNUM, bmax_initial, TXD_N など）
%   coeffs : struct（FilterCoeff, ChannelCoeff, first_C, first_F, NormRef）
%   x_sym  : 送信 PAM 列（未等化 a_k） 1×K
% 出力:
%   TX.ak, TX.bk, TX.dk など（AWGへ渡すのは通常 bk）
%
% 注意:
%   bmax は旧来どおり  bmax = bmax_initial * MODNUM / 2 で決定

K    = numel(x_sym);
M = param.MODNUM;

modtype = "PAM";
if isfield(param, "MODTYPE") && ~isempty(param.MODTYPE)
    modtype = upper(string(param.MODTYPE));
end

switch modtype
    case "QAM"
        L = sqrt(M);
        if mod(L,1) ~= 0
            error('thp_prepare_tx: MODNUM must be a perfect square for QAM.');
        end
        bmax = param.bmax_initial * (L/2);
    otherwise
        bmax = param.bmax_initial * M / 2;
end


% 入力を IIRoutput に通す（FilterCoeff/ChannelCoeff は"推定済み/橋渡し済み"を使用）
[ak,bk,ck,dk,bkfromD,bk_wowrap,ck_wowrap,ck_woTHP] = ...
    IIRoutput(K, x_sym, coeffs.FilterCoeff, coeffs.ChannelCoeff, bmax);

TX = struct();
TX.ak        = ak;          % 未等化 PAM（確認用）
TX.bk        = bk;          % ★ これを AWG へ送出（THP 後）
TX.dk        = dk;          % THP 内部系列（ログ用）
TX.bk_raw    = bk_wowrap;   % wrap 前（デバッグ用）
TX.ck_woTHP  = ck_woTHP;    % 参照用
TX.M         = M;
TX.bmax      = bmax;
end
