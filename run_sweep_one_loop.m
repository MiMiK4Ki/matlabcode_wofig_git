%% run_sweep_one_loop.m  — 極力シンプル版（reshape無し）
% 前提: simulate_case.m / get_channel_synthetic.m は既存のものを使用

%% (1) パラメータ定義（同時に param にも入れる：見た目と編集性を優先）
cutoff_coeff  = 1;                      param.cutoff_coeff  = cutoff_coeff;
SNRdB         = [0 5 10 15 20 25 30];   param.SNRdB         = SNRdB;
first_F       = 0;                      param.first_F       = first_F;
bmax_initial  = [6];                    param.bmax_initial  = bmax_initial;
MODNUM        = 2.^(2:1:3);             param.MODNUM        = MODNUM;   % [4 8]
TXD_N         = 20000;                  param.TXD_N         = TXD_N;

hFLpS         = 50;                     param.hFLpS         = hFLpS;
NoSpS         = 20;                     param.NoSpS         = NoSpS;
upconvert     = 10;                     param.upconvert     = upconvert; %#ok<NASGU>
fc            = 289.44e9;               param.fc            = fc;

Tc_but        = 1/2.16e9;               param.Tc_but        = Tc_but;
f_butterworth_cutoff = 1/Tc_but/2;      param.f_butterworth_cutoff = f_butterworth_cutoff;
calc_N        = 1000;                   param.calc_N        = calc_N;

%% (2) スイープグリッド（旧mainと同じ順序）
[cutL, snrL, fFL, bmaxL, modL] = ndgrid(cutoff_coeff, SNRdB, first_F, bmax_initial, MODNUM);
gridbox = [ numel(cutoff_coeff), numel(SNRdB), numel(first_F), numel(bmax_initial), numel(MODNUM) ];
meshN   = numel(cutL);
param.gridbox = gridbox;
param.meshloop_maxindex = meshN;

% cutoff のインデックス（チャネルの取り出しに使用）
[cutI, ~, ~, ~, ~] = ndgrid( 1:numel(cutoff_coeff), 1:numel(SNRdB), 1:numel(first_F), ...
                              1:numel(bmax_initial), 1:numel(MODNUM) ); %#ok<ASGLU>

%% (3) チャネルを cutoff ごとに1回だけ生成（使い回し）
cfg = struct('hFLpS',hFLpS,'NoSpS',NoSpS,'Tc_but',Tc_but,'fc',fc);
ChanFrames = cell(1, numel(cutoff_coeff));
for ic = 1:numel(cutoff_coeff)
    ChanFrames{ic} = get_channel_synthetic(cfg, cutoff_coeff(ic));
end

%% 推定コード

%% (4) 1本ループ（for/parfor はフラグ切替のみ）— 結果は results{i} に格納

tap = @(varargin) [];       % 研究時に中間保存しないなら no-op
results = cell(1, meshN);
rng(1);                     % 再現性を揃えたいときだけ（不要なら消してOK）

    parfor i = 1:meshN
        % この i だけの params を作って simulate_case を実行
        params = struct( ...
            'TXD_N', TXD_N, ...
            'MODNUM', modL(i), ...
            'bmax_initial', bmaxL(i), ...
            'first_F', fFL(i), ...
            'cutoff_coeff', cutL(i), ...
            'SNRdB', snrL(i), ...
            'hFLpS', hFLpS, ...
            'calc_N', calc_N, ...
            'fc', fc );
        
        results{i}= simulate_case(params, ChanFrames{cutI(i)}, tap);
        fprintf('mesh=%d / %d',i, meshN);
    end

%% reshape

R = [results{:}];
% ===== 数値（スカラー）フィールド：1行=1フィールド（超短い） =====
res.P_wTHP         = reshape([R.P_wTHP],         gridbox);
res.P_woTHP        = reshape([R.P_woTHP],        gridbox);
res.ratio          = reshape([R.ratio],          gridbox);
res.P_expand       = reshape([R.P_expand],       gridbox);
res.SINR_wTHP      = reshape([R.SINR_wTHP],      gridbox);
res.SINR_woTHP     = reshape([R.SINR_woTHP],     gridbox);
res.Capacity_wTHP  = reshape([R.Capacity_wTHP],  gridbox);
res.Capacity_woTHP = reshape([R.Capacity_woTHP], gridbox);
res.BER            = reshape([R.BER],            gridbox);
res.BER_woTHP      = reshape([R.BER_woTHP],      gridbox);
res.VARB           = reshape([R.VARB],           gridbox);
res.impulse_energ  = reshape([R.impulse_energ],  gridbox);
res.Normalized_Coeff = reshape([R.Normalized_Coeff], gridbox);
res.Hard_capacity  = reshape([R.Hard_capacity],  gridbox);

% ===== セル系（各ケースが配列/行列を持つ）：cellのまま5Dへ =====
res.r_ik         = reshape({R.r_ik},         gridbox);
% res.fitResults   = reshape({R.fitResults},   gridbox);
res.ChannelCoeff = reshape({R.ChannelCoeff}, gridbox);

% ===== 構造体（vardsi）：構造体配列を5Dへ（GUIは res.vardsi_list を参照） =====
res.vardsi_list  = reshape([R.vardsi], gridbox);

% ===== 保存（GUI 既定名に寄せる） =====
outfile = "SimResult" + num2str(param.TXD_N) + "TXDN.mat";
save(outfile, 'param', 'res', '-v7.3');
fprintf('Packed to %s\n', outfile);

