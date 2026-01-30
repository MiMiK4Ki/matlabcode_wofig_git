function IO = io_make_theory(param, ChanFrame, x_sym)
% 等化の有無に関わらず、与えられた「シンボル列 x_sym」を
% 連続系の合成インパルス応答 h_rx に通して受信波形 y_wf を生成
% 必須: ChanFrame.h_rx（*dt 済み）, ChanFrame.dt
% 入出力:
%   x_sym : 送信シンボル列（PAM; 未等化なら PAM、等化ありなら THP 後の列など）
%   IO    : .x_sym, .y_wf, .fs, .Tc, .NoSpS, .pulse(任意)

    K     = numel(x_sym);
    NoSpS = param.NoSpS;
    Tc    = param.Tc_but / param.cutoff_coeff;
    dt    = ChanFrame.dt;
    fs    = 1/dt;

    % 整合性チェック（ズレると後続の LS が崩れる）
    if abs(dt - Tc/NoSpS) > max(1e-12, 1e-6*Tc)
        warning('io_make_theory: dt(%.3g) と Tc/NoSpS(%.3g) が不一致です。param と ChanFrame の元を確認してください。', dt, Tc/NoSpS);
    end

    % 零詰めインパルス列（NoSpS ごとに x_sym を打つ）
    imp = zeros(1, (K-1)*NoSpS + 1);
    imp(1:NoSpS:end) = x_sym;

    % 連続受信波形：h_rx は *dt 済み → ここでは掛け直し不要
    y_wf = conv(imp, ChanFrame.h_rx);

    % 返却 IO
    IO = struct();
    IO.x_sym = x_sym;
    IO.y_wf  = y_wf;
    IO.fs    = fs;
    IO.Tc    = Tc;
    IO.NoSpS = NoSpS;
    if isfield(ChanFrame,'pulse'), IO.pulse = ChanFrame.pulse; else, IO.pulse = []; end
end
