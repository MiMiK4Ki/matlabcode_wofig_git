%% ================= meas_est_and_compare.m =================
clear; close all;

S = load('IQWAVEFORM_64QAM.mat');   % <- 実験チームが保存した mat
%%
% --- ここは「ファイル内の変数名」に合わせてください ---
% 送信シンボル列（QAMシンボル）
tx_sym = S.tx_upd;          % 例: tx_upd がある前提（無ければ S.tx などに変更）

% 受信IQ波形（連続サンプル列）
y_rx  = S.xout1;            % 例: xout1 がある前提（無ければ S.xout 等に変更）

% RC Tx filter パラメータ（元コードと同じはず）
beta       = 0.35;          % S.beta があれば S.beta にしてOK
nsymb      = 80;            % S.nsymb があれば S.nsymb にしてOK
filtergain = 0.3578;        % S.filtergain があれば S.filtergain にしてOK
symR       = S.symR;        % シンボルレート（無ければ手で 1e9 など）
resampling_fac = S.resampling_fac;
tx_sym = tx_sym(:);
y_rx   = y_rx(:);
NoSpS = S.nSamps;
Fs2 = symR * NoSpS * resampling_fac;    % ≒128e9
dt2 = 1 / Fs2;


%% (B) 参照波形 x_ref を Tx RC filter で生成（送信機で"既知"な波形）
rcFilt = comm.RaisedCosineTransmitFilter( ...
    'RolloffFactor', beta, ...
    'FilterSpanInSymbols', nsymb, ...
    'OutputSamplesPerSymbol', NoSpS, ...
    'Gain', filtergain);

reset(rcFilt);
x_ref = step(rcFilt, tx_sym);     % 参照Tx波形（連続サンプル列）
x_ref = S.xin1;
%% (C) 参照側の「記号中心」c_idx_ref（群遅延）
% 群遅延(サンプル) = (nsymb/2)*NoSpS
c_idx_ref = (nsymb/2) * NoSpS + 1;

% 念のため：インパルス応答のピーク位置でも確認できる
reset(rcFilt);
h_rc = step(rcFilt, [1; zeros(nsymb,1)]);   % インパルス応答（連続サンプル列）
[~, c_idx_peak] = max(abs(h_rc));
fprintf('c_idx_ref(formula)=%d, c_idx_peak(impulse)=%d\n', c_idx_ref, c_idx_peak);

%% (D) xcorr による遅延推定 + 位相推定
[xc, lags] = xcorr(y_rx, x_ref);        % 受信 vs 参照
[~, imx]   = max(abs(xc));
lag        = lags(imx);
phi        = angle(xc(imx));            % 参照に対する受信の定数位相差の推定

c_index = c_idx_ref + lag;

fprintf('lag=%d samples, phi=%.3f rad, c_index=%d\n', lag, phi, c_index);

% 位相補償：
% y_rx が (だいたい) x_ref * exp(j*theta) なら、phi ≈ theta
% => 受信を参照に合わせるなら exp(-j*phi)
y_rx_corr = y_rx * exp(-1j*phi);

% 逆に「参照側を受信に合わせる」なら
% x_ref_corr = x_ref * exp(+1j*phi);

%% (E) IO を作って LS 推定
IO = struct();
IO.y_wf   = y_rx_corr.';   % row
IO.NoSpS  = NoSpS;
IO.c_index = c_index;
IO.Tc = 1/symR;
IO.fs = symR*NoSpS;

win = struct('first_C',-2,'C_max',15,'first_F',0,'F_max',15);

coeffs = estimate_channel_ls(IO, tx_sym.', win);   % 既存の簡素版 estimate_channel_ls を使用

%% (F) 比較プロット：RC Tx filter の連続インパルス応答 vs 推定した離散タップ
% RCインパルス応答を、主タップ(=c_idx_ref)で正規化
h_rc_n = h_rc / h_rc(c_idx_ref);

% 連続線の横軸を「記号インデックス相当」に変換
m_frac = ((0:numel(h_rc_n)-1) - (c_idx_ref-1)) / NoSpS;

m_int  = win.first_C:win.C_max;

figure; hold on;
plot(m_frac, h_rc_n, 'DisplayName','|RC Tx impulse| (normalized)');
stem(m_int, coeffs.ChannelCoeff, 'filled', 'DisplayName','|Estimated ChannelCoeff|');
grid on; xlabel('symbol index m'); ylabel('magnitude');
title('RC Tx impulse (continuous) vs estimated symbol-spaced channel (discrete)');
legend('show','Location','best');

%% (G) 追加の確認：相関の形（ピークが鋭いか）
figure;
plot(lags, abs(xc)); grid on;
xlabel('lag [samples]'); ylabel('|xcorr|');
title('xcorr(y\_rx, x\_ref) magnitude');

%% (H) 追加の確認:xin1 xout1(位相補正） symbolの重ね合わせ
Kp = 300;
idx = c_index + (0:Kp-1)*NoSpS;
idx = idx(idx>=1 & idx<=numel(IO.y_wf));
y_smp = IO.y_wf(idx);


%% ================= end =================
