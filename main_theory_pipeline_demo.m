%% ================== main_theory_est_route_stubbed.m ==================
% 推定ありルートの形を保ちながら、現状は理論値で橋渡しするテンプレ
% 依存: get_channel_synthetic, io_make_theory, chan_cont_to_discrete,
%       thp_prepare_tx, evaluate_from_io, generate_PAM, 既存ユーティリティ
% （将来ON）: estimate_c_index, estimate_channel_ls

%% (0) パラメータ（ここだけ編集）
clear
param.cutoff_coeff  = 2;
param.SNRdB         = 20;
param.first_F       = 0;
param.bmax_initial  = 6;
param.MODNUM        = 8;           % 8-PAM (=2M)
param.TXD_N         = 20000;

param.hFLpS         = 50;
param.NoSpS         = 20;
param.upconvert     = 10;          %#ok<NASGU>
param.fc            = 289.44e9;

param.Tc_but        = 1/2.16e9;
param.f_butterworth_cutoff = 1/param.Tc_but/2;
param.calc_N        = 1000;

win = struct('first_C',-2,'C_max',15,'first_F',param.first_F,'F_max',15);

rng(1);  % 再現性が欲しい場合のみ

%% (1) 理論チャネル生成（A）
ChanFrame = get_channel_synthetic(param); % h_rx(*dt), dt, (pulse 任意)
figure
plot((0:numel(ChanFrame.h_rx)-1).*ChanFrame.dt,ChanFrame.h_rx)

%% (2) 訓練 PAM（C）
x_train = generate_PAM(param.MODNUM, param.TXD_N);

%% (3) 理論 I/O（等化なし）（E）
IO_wo = io_make_theory(param, ChanFrame, x_train);

%% (4) ★サンプリング中心 c_index を"推定フェーズ前"に確定
% いまは理論値を使用（主タップ最大）。将来は下の estimate_c_index をON。
[~, c_idx_theory] = max(abs(ChanFrame.h_rx));    % ← 現在はこれを採用
IO_wo.c_index = c_idx_theory;

% （将来ON）位相・記号ずれを走査して x と y だけで c_index を求める
% [c_idx_est, syncInfo] = estimate_c_index(IO_wo, x_train, 'Q',15,'Smax',300,'skipSyms',50);
% IO_wo.c_index = c_idx_est;

fprintf('[Timing] c_index=%d\n', IO_wo.c_index);

[~,c_idx_sRRC]=max(abs(ChanFrame.pulse));

ChanFrame_sRRC = ChanFrame;
ChanFrame_sRRC.h_rx = ChanFrame.pulse;

IO_woEq_sRRC = io_make_theory(param, ChanFrame_sRRC, x_train);

[xc, lags] = xcorr(IO_wo.y_wf,IO_woEq_sRRC.y_wf);

figure
plot(lags,(abs(xc)))

[~, idx] = max(abs(xc));
lag  = lags(idx);
phi1 = angle(xc(idx)); 

c_idx = c_idx_sRRC + lag;

figure; plot(xcorr(xcorr(ChanFrame.h_rx,ChanFrame.h_rx),ChanFrame.h_but)) 
%これでもc_idx_sRRCが出せる
%% 雑音付加

K = numel(x_train);

% 連続波形がちゃんとある（理論 io_make_theory / 実験波形）なら continuous
[P_wo, infoP] = io_calc_power(IO_wo, K, param, "continuous");
IO_wo = io_add_awgn(IO_wo, P_wo, param);


%% (5) ★チャネル係数（LS 推定に差し替え）
coeffs_th = chan_cont_to_discrete(param, ChanFrame, win);

coeffs_ls = estimate_channel_ls(IO_wo, x_train, win);

% 理論ブリッジと比較したい場合だけ
coeffs_th = chan_cont_to_discrete(param, ChanFrame, win);
figure; hold on;
stem(win.first_C:win.C_max, real(coeffs_th.ChannelCoeff), 'DisplayName','theory');
stem(win.first_C:win.C_max, real(coeffs_ls.ChannelCoeff), 'DisplayName','LS');
grid on; legend show; title('ChannelCoeff compare');


% 理論ブリッジと比較したいときだけ（確認用）
% coeffs_th = chan_cont_to_discrete(param, ChanFrame, win);
% figure; hold on;
% stem(win.first_C:win.C_max, real(coeffs_th.ChannelCoeff), 'DisplayName','theory (bridge)');
% stem(win.first_C:win.C_max, real(coeffs.ChannelCoeff),   'DisplayName','LS estimate');
% grid on; legend show; title('ChannelCoeff compare (real)');




%% (6) 評価（未等化, H / woTHP）
% evaluate_from_io は IO.c_index 必須。ChanFrame 引数は不要。
Res_wo = evaluate_from_io(IO_wo, coeffs, param, x_train, 'woTHP');
fprintf('[woTHP] BER=%.3e, HardCap=%.3f SoftCap=%.3f,(shift=%d, c_index=%d)\n', ...
        Res_wo.BER, Res_wo.Hard_capacity,Res_wo.Soft_capacity, Res_wo.shift_sym, Res_wo.c_index);

%% (7) THP 送信列（G）— 実験なら TX.bk を AWG に送出
TX = thp_prepare_tx(param, coeffs, x_train);

%% (8) 理論 I/O（等化後）（E）
IO_w = io_make_theory(param, ChanFrame, TX.bk);

% 等化後のサンプリング中心。まずは未等化と同じ位相を共用。
IO_w.c_index = IO_wo.c_index;

% （将来ON）等化後、改めて c_index を再推定したい場合
% [c_idx_w, infoW] = estimate_c_index(IO_w, x_train, 'Q',15,'Smax',300,'skipSyms',50);
% IO_w.c_index = c_idx_w;
%% 雑音付加

K = numel(TX.bk);

[P_w, infoPw] = io_calc_power(IO_w, K, param, "continuous");
IO_w = io_add_awgn(IO_w, P_w, param);


%% (9) 評価（等化, H / wTHP）
Res_w = evaluate_from_io(IO_w, coeffs, param, x_train, 'wTHP');
fprintf('[wTHP ] BER=%.3e, HardCap=%.3f SoftCap=%.3f(shift=%d, c_index=%d)\n', ...
        Res_w.BER, Res_w.Hard_capacity, Res_w.Soft_capacity,Res_w.shift_sym, Res_w.c_index);

%% （任意）可視化 — 連続波形と中心サンプルの重ね合わせ（先頭200記号）
% plot_io_overlay(IO_wo, ChanFrame, 200);
% plot_io_overlay(IO_w,  ChanFrame, 200);

%% ================== end of script ==================
