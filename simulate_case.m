function R = simulate_case(params, ChanFrame, tap)
% SIMULATE_CASE  1ケースの理論評価（ZF-THP）を実行
% 前提:
%   - フィルタ係数は既知（ChanFrame.FilterCoeff / ChannelCoeff / first_C / first_F / NormRef）
%   - ただし無い場合は、ChanFrame.h_rx（*dt 済み）から自動ブリッジ（chan_cont_to_discrete を呼ぶ）
%
% 入出力:
%   params: struct（旧main相当。主要: TXD_N, MODNUM, bmax_initial, first_F, cutoff_coeff, SNRdB, hFLpS, NoSpS, Tc_but, calc_N, fc）
%   ChanFrame: struct（少なくとも h_rx, dt。係数があればそれを優先）
%   tap: デバッグ用 no-op 関数ハンドル（例: @(varargin)[]）
%
% 返り値 R: 構造体（旧mainで reshape していた項目を網羅）

if nargin < 3 || isempty(tap), tap = @(varargin)[]; end

% =========================
% (0) パラメータ取り出し
% =========================
TXD_N     = params.TXD_N;
M_mod     = params.MODNUM;          % = 2M
bmax_init = params.bmax_initial;
first_F   = params.first_F;
cutoff    = params.cutoff_coeff;
SNRdB     = params.SNRdB;
hFLpS     = params.hFLpS;
NoSpS     = params.NoSpS;
Tc_but    = params.Tc_but;
calc_N    = params.calc_N;
fc        = params.fc; %#ok<NASGU>  % 今は未使用（将来拡張用）

% 時間刻みなど
Tc    = Tc_but / cutoff;
dt    = Tc / NoSpS;
fs    = 1/dt; %#ok<NASGU>

% =========================
% (1) 係数の確定（既知が無ければ自動ブリッジ）
% =========================
hasCoeffs = isfield(ChanFrame,'FilterCoeff') && isfield(ChanFrame,'ChannelCoeff') && ...
            isfield(ChanFrame,'first_C')     && isfield(ChanFrame,'first_F')     && ...
            isfield(ChanFrame,'NormRef');

if ~hasCoeffs
    % ※理論デバッグの便宜: 連続IRから離散タップにブリッジ
    ChanFrame = chan_cont_to_discrete(params, ChanFrame);
end

FilterCoeff          = ChanFrame.FilterCoeff;
ChannelCoeff         = ChanFrame.ChannelCoeff;
first_C              = ChanFrame.first_C;
first_F_value        = ChanFrame.first_F;
Normalized_Coeff_temp= ChanFrame.NormRef;
Prefilength          = numel(FilterCoeff);

tap('coeff.Filter', FilterCoeff);
tap('coeff.Channel', ChannelCoeff);

% =========================
% (2) 成形パルス（RRC）/ Butterworth（理論合成用）
% =========================
t_rrc  = -hFLpS*Tc : dt : hFLpS*Tc;
h_rrc  = srrc_filter_wosym(t_rrc, 1, Tc, 1); % r=1, B=1（旧main踏襲）
f_cut  = 1/Tc_but/2;
[h_but, ~] = generate_butterworth_filter(f_cut, hFLpS, Tc, dt);

% 合成IR（理論確認用）
impulseres_Rx = conv(h_but, h_rrc) * dt;
tap('impulse.Rx', impulseres_Rx);

% =========================
% (3) 入力（PAM）生成
% =========================
PAM_signal = generate_PAM(M_mod, TXD_N);
tap('PAM',PAM_signal);

% =========================
% (4) THP系 IIR 出力
% =========================
bmax_value = bmax_init * (M_mod/2);
[ak, bk, ck, dk, bkfromD, bk_wowrap, ck_wowrap, ck_woTHP] = ...
    IIRoutput(TXD_N, PAM_signal, FilterCoeff, ChannelCoeff, bmax_value);
tap('bk',bk);

% Mk（folded Gaussian 用）
Mk = (ak - dk) / (2*bmax_value);

% =========================
% (5) 連続畳込み（h_array）と P_expand のための項
% =========================
L0 = length(impulseres_Rx);
maxShift = (calc_N - 1) * NoSpS;
L  = L0 + maxShift;

h_array = cell(calc_N, 1);
totalWave = zeros(1, L);

bkShort = bk(1:calc_N);
for kSym = 1:calc_N
    shift_samples = (kSym - 1) * NoSpS;
    shifted = zeros(1, L);
    shifted(1+shift_samples : shift_samples+L0) = impulseres_Rx;
    h_array{kSym} = shifted;
    totalWave = totalWave + shifted .* bkShort(kSym);
end
sumcross = sum(totalWave.^2) * dt; %#ok<NASGU>

[bh_square, cross_blbk] = calc_bh_crossterm(bkShort, h_array, dt, calc_N);
h_array = []; %#ok<NASGU>

% =========================
% (6) 連続出力と電力
% =========================
[y3_wTHP,  P_wTHP_value]  = calculate_Power(bk, TXD_N, h_but, hFLpS, NoSpS, dt, h_rrc, Tc, fc);
[y3_woTHP, P_woTHP_value] = calculate_Power(ak, TXD_N, h_but, hFLpS, NoSpS, dt, h_rrc, Tc, fc);

impulse_energ = sum(impulseres_Rx.^2) * dt;
NormEnerg_wTHP  = P_wTHP_value  * ((TXD_N-1)*Tc) / impulse_energ;
NormEnerg_woTHP = P_woTHP_value * ((TXD_N-1)*Tc) / impulse_energ;
VARB = var(bk);

P_expand = sum(diag(bh_square) + cross_blbk, 'all') / ((calc_N-1)*Tc);
ratio    = P_wTHP_value / P_expand;

% =========================
% (7) 干渉ベクトルと DSI
% =========================
plmicoeff = Normalized_Coeff_temp / abs(Normalized_Coeff_temp);

% 0ラグの位置（m=0）のインデックス
idx_m0 = 1 - first_C;  % ChannelCoeff は m=first_C:... 並び
interference_woTHP = plmicoeff .* ck_woTHP(abs(first_C)+1:end) ...
                   - abs(ChannelCoeff(idx_m0)) .* ak(1:end-abs(first_C));

interference_wTHP  = ck(abs(first_C - first_F_value)+1:end) ...
                   - dk(1:end-abs(first_C - first_F_value));

epsilon    = 1.0e-5;
dk_woTHP   = ak;
ik_woTHP   = interference_woTHP;  % ← 既知の修正（1.2）
sk         = ck;
sk_woTHP   = ck_woTHP;

dsi_normal          = calculate_dsi(dk,      sk,      interference_wTHP,  P_wTHP_value,  epsilon);
dsi_norm            = calculate_dsi(dk,      sk,      interference_wTHP,  P_wTHP_value,  epsilon, Normalized_Coeff_temp);
dsi_norm_woTHP      = calculate_dsi(dk_woTHP,sk_woTHP,ik_woTHP,           P_woTHP_value, epsilon, Normalized_Coeff_temp);

% =========================
% (8) 受信推定 / SINR
% =========================
[interferencePow_wTHP, a_hat_wTHP, ck_wnoise_wTHP, SINR_wTHP] = ...
    generate_ahat(ck,       interference_wTHP, cutoff, SNRdB, Normalized_Coeff_temp, P_wTHP_value,  dsi_norm,      bmax_value);

[interferencePow_woTHP, ~, ck_wnoise_woTHP,   SINR_woTHP] = ...
    generate_ahat(ck_woTHP, interference_woTHP, cutoff, SNRdB, Normalized_Coeff_temp, P_woTHP_value, dsi_norm_woTHP, bmax_value);

% =========================
% (9) 容量 / 折り返しガウス / BER / Hard capacity
% =========================
% 容量（保存先フォルダは空文字で回避）
foldername = '';
[Capacity_wTHP, a_hat_data] = calculate_capacity(foldername, a_hat_wTHP, PAM_signal, SNRdB, first_F_value, first_C, "wTHP",  bmax_value);
[Capacity_woTHP, ~]         = calculate_capacity(foldername, ck_wnoise_woTHP, PAM_signal, SNRdB, 0,           first_C, "woTHP", bmax_value);

% folded Gaussian（視覚化/保存はしない・fit だけ保持）
fitResults = fitFoldedGaussianPerSymbol2(ck_wnoise_wTHP, PAM_signal, first_C, first_F_value, bmax_value, Mk);

% BER（旧mainと同じずらし）
TXD          = pam_to_bits(PAM_signal, M_mod);
TXD_Rx_wTHP  = pam_to_bits(a_hat_wTHP, M_mod);
TXD_Rx_woTHP = pam_to_bits(plmicoeff .* ck_wnoise_woTHP, M_mod, ChannelCoeff(idx_m0));

shiftLen_wTHP  = abs(first_C - first_F_value) * M_mod / 2;
shiftLen_woTHP = abs(first_C)                 * M_mod / 2;

BER        = sum(xor(TXD(1:end-shiftLen_wTHP),  TXD_Rx_wTHP(shiftLen_wTHP+1:end))) / numel(TXD(1:end-shiftLen_wTHP));
BER_woTHP  = sum(xor(TXD(1:end-shiftLen_woTHP), TXD_Rx_woTHP(shiftLen_woTHP+1:end))) / numel(TXD_Rx_woTHP(shiftLen_woTHP+1:end));

% Hard capacity
[TX_sym, ~] = pam_to_symbols(PAM_signal(1:end-(first_F_value-first_C)), M_mod);
[RX_sym, ~] = pam_to_symbols(a_hat_wTHP(1+first_F_value-first_C:end),   M_mod);
[r_ik, Hard_capacity] = compute_transition_and_hardcapacity(TX_sym, RX_sym, M_mod);

% =========================
% (10) 結果構造体に格納（reshape 前提のフィールド名で）
% =========================
R = struct();
R.P_wTHP            = P_wTHP_value;
R.P_woTHP           = P_woTHP_value;
R.ratio             = ratio;
R.P_expand          = P_expand;
R.SINR_wTHP         = SINR_wTHP;
R.SINR_woTHP        = SINR_woTHP;
R.Capacity_wTHP     = Capacity_wTHP;
R.Capacity_woTHP    = Capacity_woTHP;
R.BER               = BER;
R.BER_woTHP         = BER_woTHP;
R.VARB              = VARB;
R.impulse_energ     = impulse_energ;
R.Normalized_Coeff  = Normalized_Coeff_temp;
R.Hard_capacity     = Hard_capacity;

R.vardsi = struct('normal', dsi_normal, 'normalized', dsi_norm, 'normalized_woTHP', dsi_norm_woTHP);
R.r_ik   = r_ik;
R.fitResults   = fitResults;
R.ChannelCoeff = ChannelCoeff;

% 任意：必要なら内部信号も返して解析しやすく
% R.debug = struct('ak',ak,'bk',bk,'ck',ck,'dk',dk,'y3_wTHP',y3_wTHP,'y3_woTHP',y3_woTHP);

end
