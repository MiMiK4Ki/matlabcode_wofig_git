function IO = io_add_awgn(IO, P_ref, param)
% io_add_awgn  参照電力 P_ref を使って IO.y_wf に AWGN を付加
%
% 入力:
%   IO.y_wf : 連続波形（row/colどちらでもOK）
%   P_ref   : 参照電力（io_calc_power で計算して渡す）
%   param.SNRdB, param.cutoff_coeff : SNR 設定
%
% 出力:
%   IO.y_wf   : ノイズ付加後
%   IO.noise  : 付帯情報

y = IO.y_wf(:).';

% cutoff 補正（旧 generate_ahat と同じ）
if param.cutoff_coeff >= 1
    SNR_adj_dB = 0;
else
    SNR_adj_dB = 10*log10(1/param.cutoff_coeff);
end
SNR_lin = 10.^((param.SNRdB + SNR_adj_dB)/10);

% 雑音分散（1サンプルあたり）
sigma2 = P_ref / SNR_lin;

% AWGN 生成
if isreal(y)
    n = sqrt(sigma2) * randn(size(y));
else
    n = sqrt(sigma2/2) * (randn(size(y)) + 1j*randn(size(y)));
end

IO.y_wf = y + n;

IO.noise = struct( ...
    'P_ref',      P_ref, ...
    'noise_var',  sigma2, ...
    'SNR_lin',    SNR_lin, ...
    'SNRdB_eff',  10*log10(SNR_lin));

end
