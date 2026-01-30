% 1) 正解チャネルを作る
clear
cfg = struct('hFLpS',50,'NoSpS',20,'Tc_but',1/2.16e9,'fc',289.44e9);
params = struct('TXD_N',20001,'MODNUM',8 ,'bmax_initial',6,'first_F',0, ...
                'cutoff_coeff',1,'SNRdB',10,'hFLpS',50,'calc_N',1000,'fc',289.44e9);

ChanFrame = get_channel_synthetic(cfg, params.cutoff_coeff); % 例: cutoff_coeff=1.2
save chan_synthetic.mat ChanFrame

% 2) 入力作成→（理論なら）上チャネル読込／（実験なら）装置通過→出力サンプル化
%   → 本テンプレでは理論のみの例
load chan_synthetic.mat

% 3) （必要なら）入出力からチャネル推定（将来の関数）
%   ChanFrame = estimate_channel(ioRecord, cfg);

% 4) チャネルをロードし、等化・性能評価（1ケースだけ）

[tap, getTap] = make_tap(true);

tap('sanity', 42);   
res = simulate_case(params, ChanFrame, tap);
disp(res);

S = getTap();  

if isfield(S,'PAM')
    fprintf('Captured PAM length = %d\n', numel(S.PAM));
end