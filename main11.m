% function main11
%% (1) 初期化・パス設定など
% clear all
% Base_path = fullfile('/mnt','nas','homes','Shun_ISHIHARA','sim');
   % Base_path = fullfile('\\NAS37092E','homes','Shun_ISHIHARA','sim');
   % Base_path = '/mnt/NAS37092E/Shun_ISHIHARA/sim';
   % Base_path = fullfile('\\NAS37092E','homes','Shun_ISHIHARA','sim');
   Base_path = 'sim';
%% (2) パラメータ設定
cutoff_coeff     = 1;
SNRdB            = [0 5 10 15 20 25 30];
first_F          = 0;  % filter
bmax_initial     = [6];
MODNUM           = 2.^(2:1:3);  
% cutoff_coeff     = 1.2;
% SNRdB            = 0;
% first_F          = 0;  % filter
% bmax_initial     = [2];
% MODNUM           = 2.^(2);  
TXD_N            = 20000; % 送信データ数

hFLpS            = 50;    % half Filter Length per Symbol
NoSpS            = 20;    % Number of Sample per Symbol
upconvert        = 10;
fc               = 289.44e9; % upconversion frequency

Tc_but           = 1/2.16e9;
f_butterworth_cutoff = 1/Tc_but/2;

calc_N           = 1000;  % h_array作成等に使うサンプル数

%% (3) シミュレーション用グリッド生成
[cutoff_coeff_loop, SNRdB_loop, first_F_loop, bmax_loop, modnum_loop] = ...
    ndgrid(cutoff_coeff, SNRdB, first_F, bmax_initial, MODNUM);
gridbox = [ length(cutoff_coeff), length(SNRdB), length(first_F), ...
            length(bmax_initial), length(MODNUM) ];

[cutoff_coeff_loop_index, SNRdB_loop_index, first_F_loop_index, ...
    bmax_loop_index, modnum_loop_index] = ...
    ndgrid(1:length(cutoff_coeff), 1:length(SNRdB), 1:length(first_F), 1:length(bmax_initial), 1:length(MODNUM));

meshloop_maxindex = numel(ndgrid(cutoff_coeff, SNRdB, first_F, bmax_initial, MODNUM));

%% (4) チェックポイント用フォルダの設定
% 親フォルダとして TXD_N が含まれるフォルダ（例: .../30700TXD_N）を使用
parentFolder = fullfile(Base_path, num2str(TXD_N) + "TXD_N");
if ~exist(parentFolder, 'dir')
    mkdir(parentFolder);
end
checkpointFolder = fullfile(parentFolder, 'checkpoints');
if ~exist(checkpointFolder, 'dir')
    mkdir(checkpointFolder);
end

%% (5) parforループ（各反復をチェックポイント化）
for meshloop_index = 1:meshloop_maxindex
    % 各反復に対して固有のチェックポイントファイル名を設定
    % tic
    checkpointFile = fullfile(checkpointFolder, sprintf('checkpoint_%d.mat', meshloop_index));
    
    % 既に計算済みならスキップ
    if exist(checkpointFile, 'file')
        fprintf('Iteration %d already computed, skipping.\n', meshloop_index);
    else
        %%%%% (A) ループ内で使うパラメータの取り出し %%%%%
        cutoff_coeff_value = cutoff_coeff_loop(meshloop_index);
        SNRdB_value        = SNRdB_loop(meshloop_index);
        first_F_value      = first_F_loop(meshloop_index);
        Modnum_value       = modnum_loop(meshloop_index);
        bmax_value         = bmax_loop(meshloop_index) .* Modnum_value ./ 2;
        
        cutoff_coeff_index = cutoff_coeff_loop_index(meshloop_index);
        SNRdB_index        = SNRdB_loop_index(meshloop_index);
        first_F_index      = first_F_loop_index(meshloop_index);
        bmax_index         = bmax_loop_index(meshloop_index);
        
        f_cutoff = cutoff_coeff_value * f_butterworth_cutoff;
        Tc       = Tc_but / cutoff_coeff_value;
        deltat   = Tc / NoSpS;
        filter_time = -hFLpS*Tc : deltat : hFLpS*Tc;
        fs       = 1 / deltat;
        
        %%%%% (B) RRCフィルタ生成 %%%%%
        B = 1;  r = 1;
        h = srrc_filter_wosym(filter_time, r, Tc, B);
        
        % foldername0 は各反復ごとに固有の結果保存先
        foldername0 = fullfile(Base_path, num2str(TXD_N) + "TXD_N", ...
            num2str(Modnum_value) + "Modnum" + num2str(round(cutoff_coeff_value*100)) + "Ws" + num2str(r) + "r");
        if ~exist(foldername0, 'dir')
            mkdir(foldername0);
        end
        
        %%%%% (C) バタワースフィルタ生成 %%%%%
        [h_butterworth, butterworth_t] = generate_butterworth_filter(f_butterworth_cutoff, hFLpS, Tc, deltat);
        % fi = figure;
        % plot(butterworth_t, h_butterworth);
        % xlabel('timex'); ylabel('Amplitude');
        % setInterpreterLatex();
        % saveCurrentFigure(fi, foldername0, 'ButImpulse');
        % 
        % freq_but = ((0:numel(butterworth_t)-1)./numel(butterworth_t) - 0.5)*fs;
        % freq_rrc = ((0:numel(filter_time)-1)./numel(filter_time) - 0.5)*fs;
        % fi = figure;
        % semilogy(freq_but, fftshift(abs(fft(h_butterworth)).*deltat)); hold on;
        % semilogy(freq_rrc, fftshift(abs(fft(h)).*deltat));
        % xlabel('Freq'); legend('but','rrc'); setInterpreterLatex();
        % saveCurrentFigure(fi, foldername0, 'FreqRes');
        % foldername0_toc=toc
        %%%%% (D) 変調信号生成 (PAM) %%%%%
        [PAM_signal] = generate_PAM(Modnum_value, TXD_N);
        
        %%%%% (E) インパルス応答の合成（RRC+Butterworth） %%%%%
        impulseres_Tx1 = h;
        impulseres_Tx2 = conv(h_butterworth, impulseres_Tx1).*deltat;
        t = 0:deltat:(2*hFLpS*NoSpS + numel(butterworth_t) - 1)*deltat;
        impulseres_Rx = impulseres_Tx2;
        filterLength = 40;
        [c, c_index] = max(impulseres_Rx);
        sampling_delay = c_index - hFLpS*NoSpS + 1;
        sampling_index_int = -filterLength:filterLength;
        sampling_index_int_continuous = (-filterLength*NoSpS):(filterLength*NoSpS);
        sampling_index = c_index + sampling_index_int.*NoSpS;
        Sampling_impulseres_Rx = impulseres_Rx(sampling_index);
        Normalized_Sampling_ZF = Sampling_impulseres_Rx ./ Sampling_impulseres_Rx(1+filterLength);
        sampling_index_continuous = c_index + sampling_index_int_continuous;
        % fi = figure;
        % plot(sampling_index_int_continuous./NoSpS, impulseres_Rx(sampling_index_continuous)); hold on;
        % stem(sampling_index_int, Sampling_impulseres_Rx, 'o');
        % xlim([-5 12]); xlabel('Index'); ylabel('Amplitude'); setInterpreterLatex();
        % saveCurrentFigure(fi, foldername0, 'Impulse');
        
        %%%%% (F) チャネル係数の切り出し %%%%%
        first_C = -5;
        C_max = 15; F_max = 15;
        truncatedindex_G = find(sampling_index_int >= first_C & sampling_index_int <= C_max);
        truncatedindex_F = find(sampling_index_int >= first_F_value & sampling_index_int <= F_max);
        Normalized_Coeff_temp = Sampling_impulseres_Rx(truncatedindex_F(1));
        ChannelCoeff = Sampling_impulseres_Rx(truncatedindex_G) ./ Normalized_Coeff_temp;
        FilterCoeff  = Sampling_impulseres_Rx(truncatedindex_F) ./ Normalized_Coeff_temp;
        FilterCoeff(1) = FilterCoeff(1) - 1;
        Prefilength = numel(FilterCoeff);
        % fi = figure;
        % stem(first_C:C_max, ChannelCoeff); hold on;
        % stem(first_F_value:F_max, FilterCoeff);
        % xlabel('Index'); ylabel('Amplitude'); setInterpreterLatex();
        % close(fi);
        K = TXD_N;
        plmicoeff = Normalized_Coeff_temp ./ abs(Normalized_Coeff_temp);
        
        
        %%%%% (G) THP系IIR出力計算 (bk, ck, dk等) %%%%%
        foldername = fullfile(foldername0, num2str(10*bmax_value) + "bmax");
        if ~exist(foldername, 'dir')
            mkdir(foldername);
        end
        [ak, bk, ck, dk, bkfromD, bk_wowrap, ck_wowrap, ck_woTHP] = ...
            IIRoutput(K, PAM_signal, FilterCoeff, ChannelCoeff, bmax_value);
        bkShort = bk(1:calc_N);

        %%%%% 2025/5/13 NEW Mk
        
        Mk=(ak-dk)/(2*bmax_value);    

        %%%%% (H) 連続波形をまとめた h_array の構築 %%%%%
        L0 = length(impulseres_Rx);
        maxShift = (calc_N - 1) * NoSpS;
        L = L0 + maxShift;
        h_array = cell(calc_N, 1);
        totalWave = zeros(1, L);
        for kSym = 1:calc_N
            shift_samples = (kSym - 1)*NoSpS;
            shiftedWave = zeros(1, L);
            startPos = 1 + shift_samples;
            endPos = shift_samples + L0;
            shiftedWave(startPos:endPos) = impulseres_Rx;
            h_array{kSym} = shiftedWave;
            % plot(shiftedWave .* bkShort(kSym));
            totalWave = totalWave + shiftedWave .* bkShort(kSym);
        end
        sumcross = sum(totalWave.^2) * deltat;
        
        %%%%% (I) クロスターム計算 %%%%%
        [bh_square, cross_blbk] = calc_bh_crossterm(bkShort, h_array, deltat, calc_N);
        h_array=[];
        %%%%% (J) 出力波形とパワー計算 %%%%%
        [y3_wTHP,  P_wTHP_value] = calculate_Power(bk, TXD_N, h_butterworth, hFLpS, NoSpS, deltat, h, Tc, fc);
        [y3_woTHP, P_woTHP_value] = calculate_Power(ak, TXD_N, h_butterworth, hFLpS, NoSpS, deltat, h, Tc, fc);
        impulse_energ_temp = sum(impulseres_Rx.^2) * deltat;
        NormEnerg_wTHP_temp  = P_wTHP_value * ((TXD_N-1)*Tc) / impulse_energ_temp;
        NormEnerg_woTHP_temp = P_woTHP_value * ((TXD_N-1)*Tc) / impulse_energ_temp;
        VARB_temp = var(bk);
        sumy3 = sum(y3_wTHP.^2) * deltat;
        P_expand = sum(diag(bh_square) + cross_blbk, 'all') / ((calc_N-1)*Tc);
        ratio = P_wTHP_value / P_expand;
        fprintf('mesh=%d / %d: direct=%.6g, expand=%.6g, ratio=%.4f\n',...
            meshloop_index, meshloop_maxindex, P_wTHP_value, P_expand, ratio);
        % crossterm_toc = toc
        % crossterm_toc
        %%%%% (K) 受信ノイズや干渉、SINRの計算 %%%%%
        ChannelCoeff_value = ChannelCoeff(first_F_index);
        interference_woTHP = plmicoeff .* ck_woTHP(abs(first_C)+1:end) - abs(ChannelCoeff(1-first_C)) .* ak(1:end-abs(first_C));
        interference_wTHP  = ck(abs(first_C - first_F_value)+1:end) - dk(1:end-abs(first_C - first_F_value));
        TXD_N_plot = min(200, TXD_N);
        t_truncated = 0:deltat:((2*hFLpS + TXD_N_plot -1)*NoSpS)*deltat;
        t1_truncated = c_index*deltat:NoSpS*deltat:((TXD_N_plot-1)*NoSpS + c_index)*deltat;
        % saveplot_THP_Tx_Rx(foldername, t_truncated, t1_truncated, ...
        %     y3_wTHP, y3_woTHP, ak, bk, ck, dk, ...
        %     bkfromD, bk_wowrap, ck_wowrap, ck_woTHP, ...
        %     Normalized_Coeff_temp, first_C, first_F_value, first_F_index, ...
        %     ChannelCoeff_value, P_wTHP_value, P_woTHP_value, ...
        %     interference_wTHP, interference_woTHP, Tc);
        % saveplot_THP_Tx_Rx_toc = toc
        % saveplot_THP_Tx_Rx_toc
        %%%%% (L) DSIパワー計算 %%%%%
        ik = interference_wTHP;
        sk = ck;
        epsilon = 1.0e-5;
        dk_woTHP = ak;
        ik_woTHP = interference_woTHP;
        sk_woTHP = ck_woTHP;
        dsi_normal = calculate_dsi(dk, sk, ik, P_wTHP_value, epsilon);
        save_data_dsi(foldername, first_F_value, "", dsi_normal);
        dsi_norm = calculate_dsi(dk, sk, ik, P_wTHP_value, epsilon, Normalized_Coeff_temp);
        save_data_dsi(foldername, first_F_value, "Normalized", dsi_norm);
        dsi_norm_woTHP = calculate_dsi(dk_woTHP, sk_woTHP, ik_woTHP, P_woTHP_value, epsilon, Normalized_Coeff_temp);
        save_data_dsi(foldername, first_F_value, "Normalized_woTHP", dsi_norm_woTHP);
        % 通常のDSIも再計算
        dsi_normal = calculate_dsi(dk, sk, ik, P_wTHP_value, epsilon);
        dsi_norm = calculate_dsi(dk, sk, ik, P_wTHP_value, epsilon, Normalized_Coeff_temp);
        dsi_norm_woTHP = calculate_dsi(dk_woTHP, sk_woTHP, ik_woTHP, P_woTHP_value, epsilon, Normalized_Coeff_temp);
        % DSIpower_toc = toc
        % DSIpower_toc
        %%%%% (M) シンボル推定・BER計算 %%%%%
        [interferencePow_wTHP_temp, a_hat_wTHP, ck_wnoise_wTHP, SINR_wTHP_temp] = ...
            generate_ahat(ck, interference_wTHP, cutoff_coeff_value, SNRdB_value, Normalized_Coeff_temp, P_wTHP_value, dsi_norm, bmax_value);
        [interferencePow_woTHP_temp, ~, ck_wnoise_woTHP, SINR_woTHP_temp] = ...
            generate_ahat(ck_woTHP, interference_woTHP, cutoff_coeff_value, SNRdB_value, Normalized_Coeff_temp, P_woTHP_value, dsi_norm_woTHP, bmax_value);
        % 可視化はメモリ使用量削減のため省略
        % c_cnoise_toc =toc
        %%%%% (N) 容量計算 %%%%%
        [Capacity_wTHP_temp, a_hat_data] = calculate_capacity(foldername, a_hat_wTHP, PAM_signal, SNRdB_value, first_F_value, first_C, "wTHP", bmax_value);
        [Capacity_woTHP_temp, ~] = calculate_capacity(foldername, ck_wnoise_woTHP, PAM_signal, SNRdB_value, 0, first_C, "woTHP", bmax_value);
        % Capacity_soft_toc = toc
        % Capacity_soft_toc
        %%%%% folded gaussian 作成
        fitResults = fitFoldedGaussianPerSymbol2(ck_wnoise_wTHP, PAM_signal, first_C, first_F_value, bmax_value,Mk);
        saveFoldedGaussian(foldername, a_hat_data, fitResults, first_F_value, SNRdB_value, "wTHP");
        % foldedgaussian_toc = toc
        % foldedgaussian_toc
        %%%%% (O) BER計算 %%%%%
        TXD = pam_to_bits(PAM_signal, Modnum_value);
        TXD_Rx_wTHP = pam_to_bits(a_hat_wTHP, Modnum_value);
        TXD_Rx_woTHP = pam_to_bits(plmicoeff.*ck_wnoise_woTHP, Modnum_value, ChannelCoeff(1-first_C));
        shiftLen_wTHP = abs(first_C - first_F_value)*Modnum_value/2;
        shiftLen_woTHP = abs(first_C)*Modnum_value/2;
        BER_temp = sum(xor(TXD(1:end-shiftLen_wTHP), TXD_Rx_wTHP(shiftLen_wTHP+1:end))) / numel(TXD(1:end-shiftLen_wTHP));
        BER_woTHP_temp = sum(xor(TXD(1:end-shiftLen_woTHP), TXD_Rx_woTHP(shiftLen_woTHP+1:end))) / numel(TXD_Rx_woTHP(shiftLen_woTHP+1:end));
        % fi = figure;
        % stem(ak, '-b', 'LineWidth', 2); hold on;
        % stem(a_hat_wTHP(1+first_F_value-first_C:end), ':g', 'LineWidth', 2);
        % legend('a','a hat'); xlabel('Index'); ylabel('Amplitude'); xlim([0 50]);
        % setInterpreterLatex();
        % saveCurrentFigure(fi, foldername, "aahat" + "fF" + num2str(first_F_value) + "snr" + num2str(SNRdB_value));
        % BER_toc=toc

        %%%%% SER & Hard capacity 計算

[TX_sym, ~] = pam_to_symbols(PAM_signal(1:end-(first_F_value-first_C)), Modnum_value);
[RX_sym, ~] = pam_to_symbols(a_hat_wTHP(1+first_F_value-first_C:end), Modnum_value);
[r_ik_temp, Hard_capacity_temp] = compute_transition_and_hardcapacity(TX_sym, RX_sym, Modnum_value);
% Hard_capacity_toc = toc
% Hard_capacity_toc
        %%%%% (P) この反復の結果を構造体にまとめる %%%%%
        iterationRes = struct();
        iterationRes.P_wTHP_value = P_wTHP_value;
        iterationRes.P_woTHP_value = P_woTHP_value;
        iterationRes.ratio = ratio;
        iterationRes.P_expand = P_expand;
        iterationRes.sumy3 = sumy3;
        iterationRes.sumcross = sumcross;
        iterationRes.vardsi.normal = dsi_normal;
        iterationRes.vardsi.normalized = dsi_norm;
        iterationRes.vardsi.normalized_woTHP = dsi_norm_woTHP;
        iterationRes.interferencePow_wTHP = interferencePow_wTHP_temp;
        iterationRes.SINR_wTHP = SINR_wTHP_temp;
        iterationRes.interferencePow_woTHP = interferencePow_woTHP_temp;
        iterationRes.SINR_woTHP = SINR_woTHP_temp;
        iterationRes.Capacity_wTHP = Capacity_wTHP_temp;
        iterationRes.Capacity_woTHP = Capacity_woTHP_temp;
        iterationRes.BER = BER_temp;
        iterationRes.BER_woTHP = BER_woTHP_temp;
        iterationRes.impulse_energ = impulse_energ_temp;
        iterationRes.NormEnerg_wTHP = NormEnerg_wTHP_temp;
        iterationRes.NormEnerg_woTHP = NormEnerg_woTHP_temp;
        iterationRes.VARB = VARB_temp;
        iterationRes.Normalized_Coeff = Normalized_Coeff_temp;
        iterationRes.r_ik = r_ik_temp;
        iterationRes.Hard_capacity = Hard_capacity_temp;
        iterationRes.fitResults = fitResults;
        iterationRes.ChannelCoeff = ChannelCoeff;
        % チェックポイントファイルとして保存（-v7.3は大容量データ向け）
        parSave(checkpointFile, iterationRes);
        
    end
    % iteration_end  =toc
    % iteration_end


    % close all

end % parfor

%% (6) チェックポイントファイルを読み込み、各反復結果を統合
% 各変数用の配列を初期化
P_wTHP_list     = zeros(1, meshloop_maxindex);
P_woTHP_list    = zeros(1, meshloop_maxindex);
ratio_list      = zeros(1, meshloop_maxindex);
P_expand_list   = zeros(1, meshloop_maxindex);
sumy3_list      = zeros(1, meshloop_maxindex);
sumcross_list   = zeros(1, meshloop_maxindex);
vardsi_list     = repmat(struct('normal', [], 'normalized', [], 'normalized_woTHP', []), 1, meshloop_maxindex);
interferencePow_wTHP_list = zeros(1, meshloop_maxindex);
SINR_wTHP_list            = zeros(1, meshloop_maxindex);
interferencePow_woTHP_list = zeros(1, meshloop_maxindex);
SINR_woTHP_list            = zeros(1, meshloop_maxindex);
Capacity_wTHP_list        = zeros(1, meshloop_maxindex);
Capacity_woTHP_list       = zeros(1, meshloop_maxindex);
BER_list                  = zeros(1, meshloop_maxindex);
BER_woTHP_list            = zeros(1, meshloop_maxindex);
impulse_energ_list        = zeros(1, meshloop_maxindex);
NormEnerg_wTHP_list       = zeros(1, meshloop_maxindex);
NormEnerg_woTHP_list      = zeros(1, meshloop_maxindex);
VARB_list                 = zeros(1, meshloop_maxindex);
Normalized_Coeff_list     = zeros(1, meshloop_maxindex);
r_ik_list = cell(1, meshloop_maxindex); % M*M 行列のため
Hard_capacity_list = zeros(1, meshloop_maxindex);
fitResults_list    = cell(1, meshloop_maxindex);
ChannelCoeff_list = cell(1, meshloop_maxindex);
for i = 1:meshloop_maxindex
    checkpointFile = fullfile(checkpointFolder, sprintf('checkpoint_%d.mat', i));
    if exist(checkpointFile, 'file')
        S = load(checkpointFile, 'iterationRes');
        iterationRes = S.iterationRes;
        P_wTHP_list(i)     = iterationRes.P_wTHP_value;
        P_woTHP_list(i)    = iterationRes.P_woTHP_value;
        ratio_list(i)      = iterationRes.ratio;
        P_expand_list(i)   = iterationRes.P_expand;
        sumy3_list(i)      = iterationRes.sumy3;
        sumcross_list(i)   = iterationRes.sumcross;
        vardsi_list(i)     = iterationRes.vardsi;
        interferencePow_wTHP_list(i) = iterationRes.interferencePow_wTHP;
        SINR_wTHP_list(i)            = iterationRes.SINR_wTHP;
        interferencePow_woTHP_list(i) = iterationRes.interferencePow_woTHP;
        SINR_woTHP_list(i)            = iterationRes.SINR_woTHP;
        Capacity_wTHP_list(i)        = iterationRes.Capacity_wTHP;
        Capacity_woTHP_list(i)       = iterationRes.Capacity_woTHP;
        BER_list(i)                  = iterationRes.BER;
        BER_woTHP_list(i)            = iterationRes.BER_woTHP;
        impulse_energ_list(i)        = iterationRes.impulse_energ;
        NormEnerg_wTHP_list(i)       = iterationRes.NormEnerg_wTHP;
        NormEnerg_woTHP_list(i)      = iterationRes.NormEnerg_woTHP;
        VARB_list(i)                 = iterationRes.VARB;
        Normalized_Coeff_list(i)     = iterationRes.Normalized_Coeff;
        r_ik_list{i} = iterationRes.r_ik;
        Hard_capacity_list(i) = iterationRes.Hard_capacity;
        fitResults_list{i}     = iterationRes.fitResults;     
        ChannelCoeff_list{i} = iterationRes.ChannelCoeff;
    else
        warning('Missing checkpoint file for iteration %d', i);
    end
end

%% (7) 結果のreshapeと統合
res = struct();
res.P_wTHP            = reshape(P_wTHP_list, gridbox);
res.P_woTHP           = reshape(P_woTHP_list, gridbox);
res.ratio             = reshape(ratio_list, gridbox);
res.P_expand          = reshape(P_expand_list, gridbox);
res.sumy3             = reshape(sumy3_list, gridbox);
res.sumcross          = reshape(sumcross_list, gridbox);
res.vardsi_list       = reshape(vardsi_list, gridbox);
res.interferencePow_wTHP = reshape(interferencePow_wTHP_list, gridbox);
res.SINR_wTHP            = reshape(SINR_wTHP_list, gridbox);
res.interferencePow_woTHP = reshape(interferencePow_woTHP_list, gridbox);
res.SINR_woTHP           = reshape(SINR_woTHP_list, gridbox);
res.Capacity_wTHP        = reshape(Capacity_wTHP_list, gridbox);
res.Capacity_woTHP       = reshape(Capacity_woTHP_list, gridbox);
res.BER                  = reshape(BER_list, gridbox);
res.BER_woTHP            = reshape(BER_woTHP_list, gridbox);
res.impulse_energ        = reshape(impulse_energ_list, gridbox);
res.NormEnerg_wTHP       = reshape(NormEnerg_wTHP_list, gridbox);
res.NormEnerg_woTHP      = reshape(NormEnerg_woTHP_list, gridbox);
res.VARB                 = reshape(VARB_list, gridbox);
res.Normalized_Coeff     = reshape(Normalized_Coeff_list, gridbox);
res.r_ik = reshape(r_ik_list, gridbox);
res.Hard_capacity = reshape(Hard_capacity_list, gridbox);
res.fitResults =reshape(fitResults_list,gridbox);
res.ChannelCoeff = reshape(ChannelCoeff_list,gridbox);
%% (8) parfor終了後のパラメータ構造体作成
param = struct();
param.cutoff_coeff    = cutoff_coeff;
param.SNRdB           = SNRdB;
param.first_F         = first_F;
param.bmax_initial    = bmax_initial;
param.MODNUM          = MODNUM;
param.hFLpS           = hFLpS;
param.NoSpS           = NoSpS;
param.upconvert       = upconvert;
param.fc              = fc;
param.Tc_but          = Tc_but;
param.f_butterworth_cutoff = f_butterworth_cutoff;
param.TXD_N           = TXD_N;
param.calc_N          = calc_N;
param.gridbox         = gridbox;
param.meshloop_maxindex = meshloop_maxindex;

%% (9) 結果の保存
filename = fullfile(Base_path, "SimResult" + num2str(TXD_N) + "TXDN.mat");
save(filename, 'param', 'res', '-v7.3');

fprintf('Simulation completed. Final result saved to %s\n', filename);

% end