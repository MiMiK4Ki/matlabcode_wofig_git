function Out = run_one_case_3routes_noplot(param, win, ChanFrame_th, x_train, coeffs_exp, ChanFrame_s2p, param_s2p_base)
% run_one_case_3routes_noplot
% 1ケース分を「プロット無し」で実行して、結果(スカラー中心)だけ返す。
%
% ルート:
%   (A) th  : theory channel + theory IO
%   (C) px  : exp-est coeffs -> approx ChanFrame + theory IO
%   (E) s2p : s2p device + sRRC ChanFrame + theory IO
%
% 出力 Out:
%   Out.th_*, Out.px_*, Out.s2p_* それぞれ woTHP/wTHP の
%     BER / HardCap / SoftCap / EVM / EVMref

Nsym_train = numel(x_train);

% =========================================================
% (A) Theory route: channel=theory, IO=theory
% =========================================================
IO_th_wo = io_make_theory(param, ChanFrame_th, x_train);
[~, IO_th_wo.c_index] = max(abs(ChanFrame_th.h_rx));

[P_th_wo, ~] = io_calc_power(IO_th_wo, Nsym_train, param, "continuous");
IO_th_wo = io_add_awgn(IO_th_wo, P_th_wo, param);

coeffs_th_ls = estimate_channel_ls(IO_th_wo, x_train, win);
coeffs_use   = coeffs_th_ls;

Res_th_wo = evaluate_from_io(IO_th_wo, coeffs_use, param, x_train, 'woTHP');

TX_th = thp_prepare_tx(param, coeffs_use, x_train);

IO_th_w = io_make_theory(param, ChanFrame_th, TX_th.bk);
IO_th_w.c_index = IO_th_wo.c_index;

Nsym_bk = numel(TX_th.bk);
[P_th_w, ~] = io_calc_power(IO_th_w, Nsym_bk, param, "continuous");
IO_th_w = io_add_awgn(IO_th_w, P_th_w, param);

Res_th_w = evaluate_from_io(IO_th_w, coeffs_use, param, x_train, 'wTHP');

Out = struct();

% ---- th ----
Out.th_BER_woTHP      = Res_th_wo.BER;
Out.th_HardCap_woTHP  = Res_th_wo.Hard_capacity;
Out.th_SoftCap_woTHP  = Res_th_wo.Soft_capacity;
Out.th_EVM_woTHP      = Res_th_wo.EVM;
Out.th_EVMref_woTHP   = NaN; if isfield(Res_th_wo,'EVMref'), Out.th_EVMref_woTHP = Res_th_wo.EVMref; end

Out.th_BER_wTHP       = Res_th_w.BER;
Out.th_HardCap_wTHP   = Res_th_w.Hard_capacity;
Out.th_SoftCap_wTHP   = Res_th_w.Soft_capacity;
Out.th_EVM_wTHP       = Res_th_w.EVM;
Out.th_EVMref_wTHP    = NaN; if isfield(Res_th_w,'EVMref'), Out.th_EVMref_wTHP = Res_th_w.EVMref; end

% =========================================================
% (C) PartialExp route: channel=exp-est, IO=theory
% =========================================================
if nargin < 5 || isempty(coeffs_exp)
    Out.px_BER_woTHP     = NaN;
    Out.px_HardCap_woTHP = NaN;
    Out.px_SoftCap_woTHP = NaN;
    Out.px_EVM_woTHP     = NaN;
    Out.px_EVMref_woTHP  = NaN;

    Out.px_BER_wTHP      = NaN;
    Out.px_HardCap_wTHP  = NaN;
    Out.px_SoftCap_wTHP  = NaN;
    Out.px_EVM_wTHP      = NaN;
    Out.px_EVMref_wTHP   = NaN;
else
    ChanFrame_px = chan_coeffs_to_chanframe_interp(param, coeffs_exp, win, "linear");

    IO_px_wo = io_make_theory(param, ChanFrame_px, x_train);
    IO_px_wo.c_index = 1 + (0 - win.first_C) * param.NoSpS;

    [P_px_wo, ~] = io_calc_power(IO_px_wo, Nsym_train, param, "continuous");
    IO_px_wo = io_add_awgn(IO_px_wo, P_px_wo, param);

    Res_px_wo = evaluate_from_io(IO_px_wo, coeffs_exp, param, x_train, 'woTHP');

    TX_px = thp_prepare_tx(param, coeffs_exp, x_train);

    IO_px_w = io_make_theory(param, ChanFrame_px, TX_px.bk);
    IO_px_w.c_index = IO_px_wo.c_index;

    Nsym_bk2 = numel(TX_px.bk);
    [P_px_w, ~] = io_calc_power(IO_px_w, Nsym_bk2, param, "continuous");
    IO_px_w = io_add_awgn(IO_px_w, P_px_w, param);

    Res_px_w = evaluate_from_io(IO_px_w, coeffs_exp, param, x_train, 'wTHP');

    Out.px_BER_woTHP      = Res_px_wo.BER;
    Out.px_HardCap_woTHP  = Res_px_wo.Hard_capacity;
    Out.px_SoftCap_woTHP  = Res_px_wo.Soft_capacity;
    Out.px_EVM_woTHP      = Res_px_wo.EVM;
    Out.px_EVMref_woTHP   = NaN; if isfield(Res_px_wo,'EVMref'), Out.px_EVMref_woTHP = Res_px_wo.EVMref; end

    Out.px_BER_wTHP       = Res_px_w.BER;
    Out.px_HardCap_wTHP   = Res_px_w.Hard_capacity;
    Out.px_SoftCap_wTHP   = Res_px_w.Soft_capacity;
    Out.px_EVM_wTHP       = Res_px_w.EVM;
    Out.px_EVMref_wTHP    = NaN; if isfield(Res_px_w,'EVMref'), Out.px_EVMref_wTHP = Res_px_w.EVMref; end
end

% =========================================================
% (E) S2P route: channel=s2p+sRRC, IO=theory
% =========================================================
if nargin < 6 || isempty(ChanFrame_s2p) || nargin < 7 || isempty(param_s2p_base)
    Out.s2p_BER_woTHP     = NaN;
    Out.s2p_HardCap_woTHP = NaN;
    Out.s2p_SoftCap_woTHP = NaN;
    Out.s2p_EVM_woTHP     = NaN;
    Out.s2p_EVMref_woTHP  = NaN;

    Out.s2p_BER_wTHP      = NaN;
    Out.s2p_HardCap_wTHP  = NaN;
    Out.s2p_SoftCap_wTHP  = NaN;
    Out.s2p_EVM_wTHP      = NaN;
    Out.s2p_EVMref_wTHP   = NaN;
else
    % s2p側の Tc_but をケース param に反映（SNR等はケース param のまま）
    param_s2p = param;
    if isfield(param_s2p_base,'Tc_but')
        param_s2p.Tc_but = param_s2p_base.Tc_but;
    end
    if isfield(param_s2p_base,'f_butterworth_cutoff')
        param_s2p.f_butterworth_cutoff = param_s2p_base.f_butterworth_cutoff;
    end

    % s2pは「既知の連続チャネル」なので bridge を採用（=理想チャネル知識）
    coeffs_s2p = chan_cont_to_discrete(param_s2p, ChanFrame_s2p, win);

    % woTHP
    IO_s2p_wo = io_make_theory(param_s2p, ChanFrame_s2p, x_train);
    [~, IO_s2p_wo.c_index] = max(abs(ChanFrame_s2p.h_rx));

    [P_s2p_wo, ~] = io_calc_power(IO_s2p_wo, Nsym_train, param_s2p, "continuous");
    IO_s2p_wo = io_add_awgn(IO_s2p_wo, P_s2p_wo, param_s2p);

    Res_s2p_wo = evaluate_from_io(IO_s2p_wo, coeffs_s2p, param_s2p, x_train, 'woTHP');

    % wTHP
    TX_s2p = thp_prepare_tx(param_s2p, coeffs_s2p, x_train);

    IO_s2p_w = io_make_theory(param_s2p, ChanFrame_s2p, TX_s2p.bk);
    IO_s2p_w.c_index = IO_s2p_wo.c_index;

    Nsym_bk_s2p = numel(TX_s2p.bk);
    [P_s2p_w, ~] = io_calc_power(IO_s2p_w, Nsym_bk_s2p, param_s2p, "continuous");
    IO_s2p_w = io_add_awgn(IO_s2p_w, P_s2p_w, param_s2p);

    Res_s2p_w = evaluate_from_io(IO_s2p_w, coeffs_s2p, param_s2p, x_train, 'wTHP');

    Out.s2p_BER_woTHP      = Res_s2p_wo.BER;
    Out.s2p_HardCap_woTHP  = Res_s2p_wo.Hard_capacity;
    Out.s2p_SoftCap_woTHP  = Res_s2p_wo.Soft_capacity;
    Out.s2p_EVM_woTHP      = Res_s2p_wo.EVM;
    Out.s2p_EVMref_woTHP   = NaN; if isfield(Res_s2p_wo,'EVMref'), Out.s2p_EVMref_woTHP = Res_s2p_wo.EVMref; end

    Out.s2p_BER_wTHP       = Res_s2p_w.BER;
    Out.s2p_HardCap_wTHP   = Res_s2p_w.Hard_capacity;
    Out.s2p_SoftCap_wTHP   = Res_s2p_w.Soft_capacity;
    Out.s2p_EVM_wTHP       = Res_s2p_w.EVM;
    Out.s2p_EVMref_wTHP    = NaN; if isfield(Res_s2p_w,'EVMref'), Out.s2p_EVMref_wTHP = Res_s2p_w.EVMref; end
end

end
