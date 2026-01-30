function ChanFrame = get_channel_synthetic(param)
% 入力 param から理論チャネルの合成インパルス応答を生成
% 必須: param.hFLpS, param.NoSpS, param.Tc_but, param.cutoff_coeff
% 取り決め: ChanFrame.h_rx は conv(h_but, h_rrc) * dt （*dt 掛け済み）

    % 時間グリッド
    Tc = param.Tc_but / param.cutoff_coeff;
    dt = Tc / param.NoSpS;
    hFLpS = param.hFLpS;

    % RRC（送受整形パルス）
    t_rrc = -hFLpS*Tc : dt : hFLpS*Tc;
    r=param.beta;
    B=1;
    h_rrc = srrc_filter_wosym(t_rrc, r, Tc, B);  % r=1, B=1（旧来に合わせる）

    % アナログ Butterworth
    f_cut = 1/param.Tc_but/2;
    [h_but, t_but] = generate_butterworth_filter(f_cut, hFLpS, Tc, dt);

    % 合成インパルス応答（*dt 掛け）
    h_rx = conv(h_but, h_rrc) * dt;

    % 返却
    ChanFrame = struct();
    ChanFrame.h_rx  = h_rx;     % 合成IR（*dt 済み）
    ChanFrame.dt    = dt;       % サンプル刻み
    ChanFrame.pulse = h_rrc;    % 任意（マッチド用にあると便利）
    ChanFrame.meta  = struct('t_rrc',t_rrc,'t_but',t_but);
    ChanFrame.h_but = h_but;
end
