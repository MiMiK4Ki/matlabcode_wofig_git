function [Res, D] = evaluate_from_sym(y_sym, coeffs, param, ref_sym, mode)
% evaluate_from_sym
%   離散シンボル列（受信側のサンプル列）から
%     BER / HardCap / SoftCap / EVM(v2) / EVMref(v2)
%   を計算する。
%
% 入力:
%   y_sym    : 受信サンプル列（※ evaluate_from_io 側で /NormRef 済みを想定）
%   coeffs   : ChannelCoeff, first_C, first_F, NormRef など
%   param    : MODNUM, bmax_initial, SNRdB など
%   ref_sym  : 参照送信シンボル（woTHPでもwTHPでも「元のx_train」を渡す想定）
%   mode     : "woTHP" or "wTHP"
%
% 出力 Res:
%   Res.BER, Res.Hard_capacity, Res.Soft_capacity, Res.EVM, Res.EVMref
%
% 追加出力 D:
%   D.tx_al, D.rx_soft_al（整列済み）, D.rx_soft_all

arguments
    y_sym double
    coeffs struct
    param struct
    ref_sym double
    mode (1,1) string {mustBeMember(mode,["woTHP","wTHP"])}
end

y_sym  = y_sym(:).';    % row
x_ref  = ref_sym(:).';  % row

NsymAvail = min(numel(x_ref), numel(y_sym));
y_sym = y_sym(1:NsymAvail);
x_ref = x_ref(1:NsymAvail);

Mmod       = param.MODNUM;
bmax_value = param.bmax_initial * (Mmod/2);

first_C = coeffs.first_C;
first_F = coeffs.first_F;

% m=0 の位置（ChannelCoeffは m=first_C:... の並び）
idx_m0 = 1 - first_C;
A0     = coeffs.ChannelCoeff(idx_m0);

% --- modeごとの soft系列（EVM/周波数応答に使う） ---
switch mode
    case "woTHP"
        rx_soft   = (y_sym) / A0;
        shift_sym = abs(first_C);
    case "wTHP"
        rx_soft   = wrapToM(y_sym, bmax_value);
        shift_sym = abs(first_C - first_F);
end

% --- hard decision系列（HardCap/BER用） ---
[TX_sym_all, ~] = pam_to_symbols(x_ref, Mmod);
[RX_sym_all, ~] = pam_to_symbols(rx_soft, Mmod);
% figure;stem(rx_soft);
if NsymAvail <= shift_sym
    error('evaluate_from_sym: symbols too short (NsymAvail=%d, shift=%d).', NsymAvail, shift_sym);
end

TX_sym_al = TX_sym_all(1:NsymAvail-shift_sym);
RX_sym_al = RX_sym_all(shift_sym+1:NsymAvail);

% ---- Hard capacity ----
[r_ik, C_hard] = compute_transition_and_hardcapacity(TX_sym_al, RX_sym_al, Mmod);

% ---- BER ----
TX_bits = pam_to_bits(levels_from_sym(TX_sym_al, Mmod), Mmod);
RX_bits = pam_to_bits(levels_from_sym(RX_sym_al, Mmod), Mmod);
Lb = min(numel(TX_bits), numel(RX_bits));
BER = sum(xor(TX_bits(1:Lb), RX_bits(1:Lb))) / Lb;

% ---- Soft capacity ----
Soft_capacity = NaN;
if isreal(rx_soft) && isreal(x_ref)
    foldername = '';
    switch mode
        case "woTHP"
            firstF_arg = 0;
            [Soft_capacity, ~] = calculate_capacity(foldername, rx_soft, x_ref, ...
                param.SNRdB, firstF_arg, first_C, "woTHP", bmax_value);
        case "wTHP"
            firstF_arg = first_F;
            [Soft_capacity, ~] = calculate_capacity(foldername, rx_soft, x_ref, ...
                param.SNRdB, firstF_arg, first_C, "wTHP", bmax_value);
    end
end

% ---- EVM (v2) ----
tx_evm = x_ref(1:NsymAvail-shift_sym);
rx_evm = rx_soft(shift_sym+1:NsymAvail);

% EVMref は「参照星座(64QAM)」前提なので、実数PAMなどでは NaN 扱い
wantEVMref = ~(isreal(tx_evm) && isreal(rx_evm));

EVM    = NaN;
EVMref = NaN;

if exist('CalcEVMv2_silent','file') == 2
    if wantEVMref
        [~,~,EVM,EVMref] = CalcEVMv2_silent(rx_evm(:), tx_evm(:), 64);
    else
        [~,~,EVM] = CalcEVMv2_silent(rx_evm(:), tx_evm(:), 64);
        EVMref = NaN;
    end
elseif exist('CalcEVMv2','file') == 2
    % 注意: CalcEVMv2.m に表示行が残っていると sweep/parfor 時に大量出力されます
    if wantEVMref
        [~,~,EVM,EVMref] = CalcEVMv2(rx_evm(:), tx_evm(:));
    else
        [~,~,EVM] = CalcEVMv2(rx_evm(:), tx_evm(:));
        EVMref = NaN;
    end
else
    % v2が無いなら NaN
    EVM = NaN; EVMref = NaN;
end

% ---- 返り値 ----
Res = struct();
Res.BER           = BER;
Res.Hard_capacity = C_hard;
Res.Soft_capacity = Soft_capacity;
Res.EVM           = EVM;
Res.EVMref        = EVMref;
Res.r_ik          = r_ik;
Res.shift_sym     = shift_sym;

if nargout >= 2
    D = struct();
    D.tx_al       = tx_evm(:).';
    D.rx_soft_al  = rx_evm(:).';
    D.rx_soft_all = rx_soft(:).';
    D.NsymAvail   = NsymAvail;
else
    D = [];
end
end

function lv = levels_from_sym(sym_idx, Mmod)
lv_ideal = (-(Mmod-1):2:(Mmod-1));
lv = lv_ideal(sym_idx);
end
