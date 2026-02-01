function plot_txrx_overlay(IO, coeffs, param, x_train, mode, tag, K_plot, tx_dk, lin_conv_output, K_hist, saveOpt)
% plot_txrx_overlay
% 送信シンボル・受信シンボル結果・連続受信波形（正規化）を重ね描きして、
% サンプリングタイミングとズレ(first_C/first_F)の整合を目視確認する。
%
% [ADD]
%   - saveOpt を渡すと、この関数内の図を「意味ある名前」で保存（未指定は保存OFF）
%   - saveOpt は最後に渡すのを推奨。ただし互換のため、
%       tx_dk / lin_conv_output / K_hist の位置に struct/string/logical を置いても解釈する

    if nargin < 6 || isempty(tag),    tag = ""; end
    if nargin < 7 || isempty(K_plot), K_plot = min(2000, numel(x_train)); end
    if nargin < 8, tx_dk = []; end
    if nargin < 9, lin_conv_output = []; end
    if nargin < 10, K_hist = []; end
    if nargin < 11, saveOpt = []; end

    % ---- 互換: saveOpt が別位置に置かれているケースを救済 ----
    if isempty(saveOpt)
        if txrx_is_saveopt(K_hist), saveOpt = K_hist; K_hist = []; end
        if isempty(saveOpt) && txrx_is_saveopt(lin_conv_output), saveOpt = lin_conv_output; lin_conv_output = []; end
        if isempty(saveOpt) && txrx_is_saveopt(tx_dk), saveOpt = tx_dk; tx_dk = []; end
    end
    opt = txrx_parse_saveopt(saveOpt, param, coeffs, tag, mode);

    % 保存トラッカ
    figList  = gobjects(0);
    nameList = strings(0);

    % ---- 互換: 9番目がスカラーで10番目が無い場合は K_hist とみなす ----
    if isempty(K_hist)
        if ~isempty(lin_conv_output) && isscalar(lin_conv_output)
            K_hist = lin_conv_output;
            lin_conv_output = [];
        else
            K_hist = min(20000, numel(x_train));
        end
    else
        if isempty(K_hist) || ~isfinite(K_hist) || K_hist <= 0
            K_hist = min(20000, numel(x_train));
        end
    end

    mode = string(mode);

    y       = IO.y_wf(:).';
    NoSpS   = IO.NoSpS;
    c_index = IO.c_index;

    first_C = coeffs.first_C;
    first_F = coeffs.first_F;
    NormRef = coeffs.NormRef;

    idx_m0 = 1 - first_C;                 
    A0     = coeffs.ChannelCoeff(idx_m0);

    % ===== 変調種別判定 =====
    modtype = "PAM";
    if isfield(param, "MODTYPE") && ~isempty(param.MODTYPE)
        modtype = upper(string(param.MODTYPE));
    elseif ~isreal(x_train)
        modtype = "QAM";
    end

    % ===== bmax =====
    Mmod = param.MODNUM;
    switch modtype
        case "QAM"
            L = sqrt(double(Mmod));
            if abs(L - round(L)) > 1e-12
                error('plot_txrx_overlay: MODNUM must be a perfect square for QAM.');
            end
            bmax_value = param.bmax_initial * (L/2);
        otherwise
            bmax_value = param.bmax_initial * (Mmod/2);
    end

    % ===== shift & scale =====
    if mode == "woTHP"
        shift_sym = abs(first_C);
        scale = (NormRef * A0);
    else
        shift_sym = abs(first_C - first_F);
        scale = NormRef;
    end

    yN = y ./ scale;

    Krx_max = min(numel(x_train), floor((numel(yN) - c_index) / NoSpS) + 1);

    % ============================================================
    % (1) overlay (K_plot)
    % ============================================================
    K_show = min(K_plot, numel(x_train));
    K_use  = min(Krx_max, K_show + shift_sym);

    idx_rx = c_index + (first_C + (0:K_use-1)) * NoSpS;
    y_smp  = yN(idx_rx);

    switch mode
        case "woTHP"
            rx_soft = y_smp;
        case "wTHP"
            rx_soft = wrapToM(y_smp, bmax_value);
        otherwise
            error("mode must be 'woTHP' or 'wTHP'.");
    end

    % hard decision
    switch modtype
        case "QAM"
            [RX_sym_idx, const] = qam_to_symbols(rx_soft, Mmod, x_train);
            rx_hard = const(RX_sym_idx);
        otherwise
            [RX_sym_idx, ~] = pam_to_symbols(real(rx_soft), Mmod);
            rx_hard = levels_from_sym(RX_sym_idx, Mmod);
    end

    Kc = K_use - shift_sym;
    if ~isempty(tx_dk)
        tx_dk = tx_dk(:).';
        Kc = min(Kc, numel(tx_dk));
    end

    idx_common      = idx_rx(shift_sym+1 : shift_sym+Kc);
    tx_train_common = x_train(1:Kc);

    y_smp_common    = y_smp(shift_sym+1 : shift_sym+Kc);
    rx_soft_common  = rx_soft(shift_sym+1 : shift_sym+Kc);
    rx_hard_common  = rx_hard(shift_sym+1 : shift_sym+Kc);

    if ~isempty(tx_dk), tx_dk_common = tx_dk(1:Kc); else, tx_dk_common = []; end

    n1 = max(1, idx_rx(1) - 2*NoSpS);
    n2 = min(numel(yN), idx_rx(end) + 2*NoSpS);
    n  = n1:n2;

    figH = local_plot_component("Real", @real, n, yN, idx_common, y_smp_common, rx_soft_common, rx_hard_common, ...
        tx_train_common, tx_dk_common, lin_conv_output, tag, mode, Kc, shift_sym);
    [figList, nameList] = txrx_push(figList, nameList, figH, "OVERLAY_REAL");

    if ~isreal(yN) || ~isreal(tx_train_common) || ~isreal(rx_soft_common)
        figH = local_plot_component("Imag", @imag, n, yN, idx_common, y_smp_common, rx_soft_common, rx_hard_common, ...
            tx_train_common, tx_dk_common, lin_conv_output, tag, mode, Kc, shift_sym);
        [figList, nameList] = txrx_push(figList, nameList, figH, "OVERLAY_IMAG");
    end

    % ============================================================
    % (2) histogram (K_hist)
    % ============================================================
    K_hist = min(round(K_hist), numel(x_train));
    K_use_h = min(Krx_max, K_hist + shift_sym);

    if K_use_h <= shift_sym
        warning('plot_txrx_overlay: K_hist too small after shift. Skip histogram/constellation.');
        txrx_save_figlist(figList, nameList, opt, tag, mode, "plot_txrx_overlay");
        return;
    end

    idx_rx_h = c_index + (first_C + (0:K_use_h-1)) * NoSpS;
    y_smp_h  = yN(idx_rx_h);

    switch mode
        case "woTHP"
            rx_soft_h = y_smp_h;
        case "wTHP"
            rx_soft_h = wrapToM(y_smp_h, bmax_value);
    end

    Kc_h = K_use_h - shift_sym;

    [figPdf, figErr] = local_plot_hist("Real", real(x_train(1:Kc_h)).', real(rx_soft_h(shift_sym+1 : shift_sym+Kc_h)).', tag, mode);
    [figList, nameList] = txrx_push(figList, nameList, figPdf, "HIST_REAL_PDF");
    [figList, nameList] = txrx_push(figList, nameList, figErr, "HIST_REAL_ERRPDF");

    if ~isreal(rx_soft_h) || ~isreal(x_train)
        [figPdf, figErr] = local_plot_hist("Imag", imag(x_train(1:Kc_h)).', imag(rx_soft_h(shift_sym+1 : shift_sym+Kc_h)).', tag, mode);
        [figList, nameList] = txrx_push(figList, nameList, figPdf, "HIST_IMAG_PDF");
        [figList, nameList] = txrx_push(figList, nameList, figErr, "HIST_IMAG_ERRPDF");
    end

    % ============================================================
    % (3) constellation
    % ============================================================
    tx_sc = x_train(1:Kc_h);
    rx_sc = rx_soft_h(shift_sym+1 : shift_sym+Kc_h);
    figH = local_plot_constellation(tx_sc, rx_sc, tag, mode, modtype, Mmod);
    [figList, nameList] = txrx_push(figList, nameList, figH, "CONSTELLATION_IQ");

    % ============================================================
    % [ADD] save
    % ============================================================
    txrx_save_figlist(figList, nameList, opt, tag, mode, "plot_txrx_overlay");
end

function lv = levels_from_sym(sym_idx, Mmod)
    lv_ideal = (-(Mmod-1):2:(Mmod-1));
    lv = lv_ideal(sym_idx);
end

function [figList, nameList] = txrx_push(figList, nameList, figH, label)
    if isempty(figH) || ~isvalid(figH), return; end
    figList(end+1)  = figH;
    nameList(end+1) = string(label);
end

function figH = local_plot_component(tag_comp, comp, n, yN, idx_common, y_smp_common, rx_soft_common, rx_hard_common, ...
    tx_train_common, tx_dk_common, lin_conv_output, tag, mode, Kc, shift_sym)
    figH = figure; hold on;

    plot(n, comp(yN(n)), 'DisplayName',sprintf('y\\_wf / scale (%s)', tag_comp));
    plot(idx_common, comp(y_smp_common), '.', 'DisplayName',sprintf('samples (y\\_smp, %s)', tag_comp));
    plot(idx_common, comp(rx_soft_common), 'x', 'DisplayName',sprintf('RX soft (%s)', tag_comp));
    stem(idx_common, comp(rx_hard_common), 'filled', 'DisplayName',sprintf('RX hard (%s)', tag_comp));
    stem(idx_common, comp(tx_train_common), 'DisplayName',sprintf('TX x\\_train (%s)', tag_comp));

    if ~isempty(tx_dk_common)
        stem(idx_common, comp(tx_dk_common), 'DisplayName',sprintf('TX dk (%s)', tag_comp));
    end

    if ~isempty(lin_conv_output)
        stem(idx_common, comp(lin_conv_output(1:numel(idx_common))), 'DisplayName',sprintf('conv(x,h) (%s)', tag_comp));
    end

    grid on;
    xlabel('sample index');
    ylabel('amplitude (normalized)');
    title(sprintf('%s  [%s]  %s (K=%d, shift=%d)', tag, mode, tag_comp, Kc, shift_sym));
    legend('show','Location','best');
end

function [figPdf, figErr] = local_plot_hist(tag_comp, tx_h, rx_h, tag, mode)
    mu = mean(rx_h);
    sg = std(rx_h);
    if ~isfinite(sg) || (sg < 1e-12)
        sg = max(abs(rx_h - mu)) + 1e-12;
    end
    rmin = mu - 6*sg;
    rmax = mu + 6*sg;
    if rmin == rmax
        rmin = mu - 1;
        rmax = mu + 1;
    end
    Nb = 80;
    edges = linspace(rmin, rmax, Nb+1);

    lv = unique(tx_h(:));
    lv = sort(lv);

    figPdf = figure; hold on; grid on;
    for ii = 1:numel(lv)
        mask = (tx_h == lv(ii));
        data = rx_h(mask);
        if isempty(data), continue; end
        histogram(data, edges, ...
            'Normalization','pdf', ...
            'DisplayStyle','stairs', ...
            'LineWidth',1.5, ...
            'DisplayName',sprintf('Tx=%g (N=%d)', lv(ii), sum(mask)));
    end

    xlabel(sprintf('RX soft amplitude (%s)', tag_comp));
    ylabel('pdf');
    title(sprintf('%s  [%s]  RX soft pdf per Tx (%s)  (N=%d)', ...
        tag, mode, tag_comp, numel(rx_h)));
    legend('show','Location','best');

    e = rx_h(:) - tx_h(:);
    figErr = figure;
    histogram(e, 'Normalization','pdf');
    grid on;
    xlabel(sprintf('RX\\_soft - TX (%s)', tag_comp));
    ylabel('pdf');
    title(sprintf('%s  [%s]  Error pdf (%s, var=%.6g)  (N=%d)', ...
        tag, mode, tag_comp, var(e), numel(e)));
end

function figH = local_plot_constellation(tx_sc, rx_sc, tag, mode, modtype, Mmod)
    tx = tx_sc(:);
    rx = rx_sc(:);

    ok = isfinite(real(tx)) & isfinite(imag(tx)) & isfinite(real(rx)) & isfinite(imag(rx));
    tx = tx(ok);
    rx = rx(ok);

    if isempty(rx)
        warning('local_plot_constellation: empty after finite check.');
        figH = [];
        return;
    end

    N = numel(rx);
    Nmax = 8000;
    if N > Nmax
        idx = unique(round(linspace(1, N, Nmax)));
    else
        idx = 1:N;
    end

    tx_u = unique(tx);

    figH = figure; hold on; grid on; axis equal;
    scatter(real(rx(idx)), imag(rx(idx)), 10, 'filled', ...
        'DisplayName', sprintf('RX soft (%d/%d)', numel(idx), N));
    plot(real(tx_u), imag(tx_u), 'x', 'LineWidth', 2, 'MarkerSize', 10, ...
        'DisplayName', sprintf('TX ref points (unique=%d)', numel(tx_u)));

    xlabel('I');
    ylabel('Q');
    title(sprintf('%s  [%s]  IQ constellation (mod=%s, M=%d)', tag, mode, modtype, Mmod));
    legend('show','Location','best');
end

%% =========================
%% [ADD] save helpers (txrx)
%% =========================
function tf = txrx_is_saveopt(x)
    tf = false;
    if isempty(x), return; end
    if islogical(x) && isscalar(x), tf = true; return; end
    if ischar(x) || (isstring(x) && isscalar(x)), tf = true; return; end
    if isstruct(x)
        tf = isfield(x,'enable') || isfield(x,'doSave') || isfield(x,'folder') || isfield(x,'outDir') || isfield(x,'formats') || isfield(x,'dpi');
    end
end

function opt = txrx_parse_saveopt(saveOptArg, param, coeffs, tag, mode)
    opt = struct();
    opt.enable = false;
    opt.folder = '';
    opt.prefix = 'txrx';
    opt.formats = {'png','fig'};
    opt.dpi = 200;

    if nargin < 1 || isempty(saveOptArg)
        return;
    end

    if islogical(saveOptArg) && isscalar(saveOptArg)
        opt.enable = logical(saveOptArg);
    elseif ischar(saveOptArg) || (isstring(saveOptArg) && isscalar(saveOptArg))
        opt.enable = true;
        opt.folder = char(saveOptArg);
    elseif isstruct(saveOptArg)
        if isfield(saveOptArg,'enable'), opt.enable = logical(saveOptArg.enable); else, opt.enable = true; end
        if isfield(saveOptArg,'doSave'), opt.enable = logical(saveOptArg.doSave); end

        if isfield(saveOptArg,'folder'), opt.folder = char(string(saveOptArg.folder)); end
        if isfield(saveOptArg,'outDir'), opt.folder = char(string(saveOptArg.outDir)); end

        if isfield(saveOptArg,'prefix'), opt.prefix = char(string(saveOptArg.prefix)); end
        if isfield(saveOptArg,'formats'), opt.formats = saveOptArg.formats; end
        if isfield(saveOptArg,'dpi'), opt.dpi = double(saveOptArg.dpi); end
    else
        return;
    end

    if ~opt.enable
        return;
    end

    if isempty(opt.folder)
        opt.folder = txrx_make_default_folder(param, coeffs, tag, mode);
    end
    if ~exist(opt.folder,'dir'), mkdir(opt.folder); end
end

function outDir = txrx_make_default_folder(param, coeffs, tag, mode)
    ts = datestr(now,'yyyymmdd_HHMMSS');
    modtype = "NA";
    if isfield(param,'MODTYPE') && ~isempty(param.MODTYPE)
        modtype = upper(string(param.MODTYPE));
    end

    tokens = [
        "txrx"
        "mod_" + modtype
        "M_" + string(param.MODNUM)
        "SNR_" + txrx_numtag(param.SNRdB)
        "NoSpS_" + string(param.NoSpS)
        "N_" + string(param.TXD_N)
        "bmax_" + txrx_numtag(param.bmax_initial)
        "fC_" + string(coeffs.first_C)
        "fF_" + string(coeffs.first_F)
        "mode_" + string(mode)
        "tag_" + string(tag)
    ];

    name = ts + "__" + strjoin(tokens','__');
    name = txrx_sanitize(name);
    outDir = fullfile(pwd, 'results_plots', name);
end

function txrx_save_figlist(figList, nameList, opt, tag, mode, funcName)
    if ~isstruct(opt) || ~isfield(opt,'enable') || ~opt.enable
        return;
    end
    if isempty(figList), return; end

    baseTag  = txrx_sanitize(string(tag));
    baseMode = txrx_sanitize(string(mode));

    for i = 1:numel(figList)
        figH = figList(i);
        if ~isvalid(figH), continue; end

        plotKind = txrx_sanitize(nameList(i));
        base = sprintf('%s__%s__%s__%s__%02d__%s', opt.prefix, funcName, baseMode, baseTag, i, plotKind);
        base = txrx_sanitize(base);
        txrx_save_one_figure(figH, opt, base);
    end
end

function txrx_save_one_figure(figH, opt, baseName)
    fmts = opt.formats;
    if isstring(fmts), fmts = cellstr(fmts); end

    for k = 1:numel(fmts)
        fmt = lower(string(fmts{k}));
        outBase = fullfile(opt.folder, baseName);
        switch fmt
            case "fig"
                try, savefig(figH, outBase + ".fig"); catch, end
            case "png"
                try
                    if exist('exportgraphics','file') == 2
                        exportgraphics(figH, outBase + ".png", 'Resolution', opt.dpi);
                    else
                        print(figH, outBase + ".png", '-dpng', sprintf('-r%d', opt.dpi));
                    end
                catch
                end
            case "pdf"
                try
                    if exist('exportgraphics','file') == 2
                        exportgraphics(figH, outBase + ".pdf");
                    else
                        print(figH, outBase + ".pdf", '-dpdf');
                    end
                catch
                end
        end
    end
end

function s = txrx_numtag(x)
    if isempty(x) || ~isfinite(x)
        s = "nan";
        return;
    end
    s = string(sprintf('%.6g', double(x)));
    s = strrep(s, "-", "m");
    s = strrep(s, "+", "");
    s = strrep(s, ".", "p");
end

function s = txrx_sanitize(s)
    s = string(s);
    s = regexprep(s, '[\\/:*?"<>|]', '_');
    s = regexprep(s, '\s+', '_');
    s = regexprep(s, '_+', '_');
    s = strip(s, '_');
end
