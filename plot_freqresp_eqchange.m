function [h_wo,h_w] = plot_freqresp_eqchange(D_wo, D_w, win, param, figTitle, coeffs_ch, saveOpt)
% plot_freqresp_eqchange
%   等化前(woTHP)と等化後(wTHP)で「有効チャネル」をLS推定し、
%   周波数応答（Magnitude + Phase）を比較表示する。
%
% [ADD]
%   - saveOpt を渡すと、生成した図を意味ある名前で保存（未指定は保存OFF）
%   - saveOpt.plotChannelCoeffTaps=true で ChannelCoeff(正規化)タップ図も追加で出す

if nargin < 7, saveOpt = []; end
opt = freq_parse_saveopt(saveOpt, param, win, figTitle);

% 保存用トラッカ（figure作成時にラベル付け）
figList  = gobjects(0);
nameList = strings(0);

mRange = win.first_C:win.C_max;

% --- 有効チャネル(woTHP / wTHP) をLS推定 ---
h_wo = estimate_fir_ls_sym(D_wo.tx_al, D_wo.rx_soft_al, mRange);
h_w  = estimate_fir_ls_sym(D_w.tx_al,  D_w.rx_soft_al,  mRange);

% m=0 で正規化して形を比較（主タップ=1）
idx0 = 1 - mRange(1);  % m=0 が入っている位置
h_wo = h_wo / h_wo(idx0);
h_w  = h_w  / h_w(idx0);

% --- 周波数軸（シンボルレート基準） ---
Tc = param.Tc_but / param.cutoff_coeff;
fs_sym = 1/Tc;

Nfft = 4096;
f = (-Nfft/2:Nfft/2-1)/Nfft * fs_sym;

H_wo = fftshift(fft(h_wo, Nfft));
H_w  = fftshift(fft(h_w,  Nfft));

% (任意) 物理チャネル
has_ch = (nargin >= 6 && ~isempty(coeffs_ch));
if has_ch
    g = coeffs_ch.ChannelCoeff(:).';
    if numel(g) ~= numel(mRange)
        has_ch = false;
    else
        g  = g / g(idx0);
        Hg = fftshift(fft(g, Nfft));
    end
end

% =========================================================
% (1) Magnitude [dB]
% =========================================================
figH = figure; hold on;
plot(f/1e9, 20*log10(abs(H_wo)+1e-12), 'DisplayName','effective woTHP (LS)');
plot(f/1e9, 20*log10(abs(H_w )+1e-12), 'DisplayName','effective wTHP  (LS)');
if has_ch
    plot(f/1e9, 20*log10(abs(Hg)+1e-12), 'DisplayName','ChannelCoeff (given)');
end
grid on;
xlabel('Frequency [GHz]'); ylabel('Magnitude [dB]');
title(['Freq response MAG (before/after eq): ' char(figTitle)]);
legend('show','Location','best');
[figList, nameList] = freq_push(figList, nameList, figH, "MAG_dB");

% =========================================================
% (2) Phase [rad]  (unwrap)
% =========================================================
figH = figure; hold on;
plot(f/1e9, unwrap(angle(H_wo)), 'DisplayName','phase woTHP (LS)');
plot(f/1e9, unwrap(angle(H_w )), 'DisplayName','phase wTHP  (LS)');
if has_ch
    plot(f/1e9, unwrap(angle(Hg)), 'DisplayName','phase ChannelCoeff (given)');
end
grid on;
xlabel('Frequency [GHz]'); ylabel('Phase [rad] (unwrap)');
title(['Freq response PHASE (before/after eq): ' char(figTitle)]);
legend('show','Location','best');
[figList, nameList] = freq_push(figList, nameList, figH, "PHASE_rad");

% =========================================================
% (3) Group delay [ps]
% =========================================================
w = 2*pi*f;  % [rad/s]
phi_wo = unwrap(angle(H_wo));
phi_w  = unwrap(angle(H_w));

tau_wo = -gradient(phi_wo, w);  % [s]
tau_w  = -gradient(phi_w , w);  % [s]

figH = figure; hold on;
plot(f/1e9, tau_wo*1e12, 'DisplayName','group delay woTHP (LS)');
plot(f/1e9, tau_w *1e12, 'DisplayName','group delay wTHP  (LS)');
if has_ch
    phi_g = unwrap(angle(Hg));
    tau_g = -gradient(phi_g, w);
    plot(f/1e9, tau_g*1e12, 'DisplayName','group delay ChannelCoeff (given)');
end
grid on;
xlabel('Frequency [GHz]'); ylabel('Group delay [ps]');
title(['Group delay (from phase slope): ' char(figTitle)]);
legend('show','Location','best');
[figList, nameList] = freq_push(figList, nameList, figH, "GROUP_DELAY_ps");

% =========================================================
% (4) time-domain taps (Real/Imag)
% =========================================================
figH = figure; hold on;
stem(mRange , real(h_wo), 'DisplayName','h\_eff woTHP (real)');
stem(mRange + param.first_F, real(h_w ), 'DisplayName','h\_eff wTHP  (real)');
grid on; xlabel('tap index m'); ylabel('tap value');
title(['Effective taps (LS) [Real]: ' char(figTitle)]);
legend('show','Location','best');
[figList, nameList] = freq_push(figList, nameList, figH, "TAPS_REAL");

if ~isreal(h_wo) || ~isreal(h_w)
    figH = figure; hold on;
    stem(mRange , imag(h_wo), 'DisplayName','h\_eff woTHP (imag)');
    stem(mRange + param.first_F, imag(h_w ), 'DisplayName','h\_eff wTHP  (imag)');
    grid on; xlabel('tap index m'); ylabel('tap value');
    title(['Effective taps (LS) [Imag]: ' char(figTitle)]);
    legend('show','Location','best');
    [figList, nameList] = freq_push(figList, nameList, figH, "TAPS_IMAG");
end

print_metrics = @(name, hvec) local_print_metrics(name, hvec, mRange);
print_metrics("h_wo (LS)", h_wo);
print_metrics("h_w  (LS)", h_w);
if has_ch, print_metrics("ChannelCoeff (given)", g); end

% =========================================================
% (5) [ADD] ChannelCoeff taps plot (normalized by m=0)
% =========================================================
doChTapPlot = freq_getflag(opt,'plotChannelCoeffTaps', false);
if doChTapPlot && ~isempty(coeffs_ch) && isfield(coeffs_ch,'ChannelCoeff') && ~isempty(coeffs_ch.ChannelCoeff)

    gplot = coeffs_ch.ChannelCoeff(:).';
    if isfield(coeffs_ch,'first_C') && ~isempty(coeffs_ch.first_C)
        mplot = coeffs_ch.first_C + (0:numel(gplot)-1);
    else
        mplot = (win.first_C) + (0:numel(gplot)-1);
    end

    idx0g = find(mplot==0, 1);
    if ~isempty(idx0g) && abs(gplot(idx0g)) > 0
        gplot_n = gplot / gplot(idx0g);
    else
        gplot_n = gplot;
    end

    figH = figure; hold on;
    stem(mplot, real(gplot_n), 'filled', 'DisplayName','ChannelCoeff real');
    grid on; xlabel('tap index m'); ylabel('value');
    title(['ChannelCoeff taps (norm @ m=0) [Real]: ' char(figTitle)]);
    legend('show','Location','best');
    [figList, nameList] = freq_push(figList, nameList, figH, "CHCOEFF_TAPS_REAL");

    if ~isreal(gplot_n)
        figH = figure; hold on;
        stem(mplot, imag(gplot_n), 'filled', 'DisplayName','ChannelCoeff imag');
        grid on; xlabel('tap index m'); ylabel('value');
        title(['ChannelCoeff taps (norm @ m=0) [Imag]: ' char(figTitle)]);
        legend('show','Location','best');
        [figList, nameList] = freq_push(figList, nameList, figH, "CHCOEFF_TAPS_IMAG");
    end
end

% =========================================================
% [ADD] save all labeled figures
% =========================================================
freq_save_figlist(figList, nameList, opt, figTitle, "plot_freqresp_eqchange");
end


function local_print_metrics(name, h, m)
    H0  = sum(h);
    Hpi = sum(h .* ((-1).^m));
    Se  = sum(h(mod(m,2)==0));
    So  = sum(h(mod(m,2)~=0));
    fprintf('[%s] H(0)=%.4f, |H(pi)|=%.4f, ratio(pi/0)=%.2f dB\n', ...
        name, real(H0), abs(Hpi), 20*log10(abs(Hpi)/abs(H0)));
    fprintf('      Se=%.4f, So=%.4f  => H0=Se+So=%.4f, Hpi=Se-So=%.4f\n', ...
        real(Se), real(So), real(Se+So), real(Se-So));
end

%% =========================
%% [ADD] save helpers (freq)
%% =========================
function [figList, nameList] = freq_push(figList, nameList, figH, label)
    figList(end+1)  = figH;
    nameList(end+1) = string(label);
end

function opt = freq_parse_saveopt(saveOptArg, param, win, figTitle)
    opt = struct();
    opt.enable = false;
    opt.folder = '';
    opt.prefix = 'freq';
    opt.formats = {'png','fig'};
    opt.dpi = 200;
    opt.plotChannelCoeffTaps = [];

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
        if isfield(saveOptArg,'plotChannelCoeffTaps'), opt.plotChannelCoeffTaps = logical(saveOptArg.plotChannelCoeffTaps); end
    else
        return;
    end

    if ~opt.enable
        return;
    end

    if isempty(opt.folder)
        opt.folder = freq_make_default_folder(param, win, figTitle);
    end
    if ~exist(opt.folder,'dir'), mkdir(opt.folder); end

    if isempty(opt.plotChannelCoeffTaps)
        opt.plotChannelCoeffTaps = true; % 保存ONなら taps 図も基本ON
    end
end

function tf = freq_getflag(opt, field, defaultVal)
    tf = defaultVal;
    if isstruct(opt) && isfield(opt, field) && ~isempty(opt.(field))
        tf = logical(opt.(field));
    end
end

function outDir = freq_make_default_folder(param, win, figTitle)
    ts = datestr(now,'yyyymmdd_HHMMSS');
    modtype = "NA";
    if isfield(param,'MODTYPE') && ~isempty(param.MODTYPE)
        modtype = upper(string(param.MODTYPE));
    end

    tokens = [
        "freqresp"
        "mod_" + modtype
        "M_"   + string(param.MODNUM)
        "SNR_" + freq_numtag(param.SNRdB)
        "NoSpS_" + string(param.NoSpS)
        "N_" + string(param.TXD_N)
        "bmax_" + freq_numtag(param.bmax_initial)
        "C_" + string(win.first_C) + "to" + string(win.C_max)
        "F_" + string(win.first_F) + "to" + string(win.F_max)
        "tag_" + string(figTitle)
    ];

    name = ts + "__" + strjoin(tokens','__');
    name = freq_sanitize(name);

    outDir = fullfile(pwd, 'results_plots', name);
end

function freq_save_figlist(figList, nameList, opt, figTitle, funcName)
    if ~isstruct(opt) || ~isfield(opt,'enable') || ~opt.enable
        return;
    end
    if isempty(figList), return; end

    baseTitle = freq_sanitize(string(figTitle));

    for i = 1:numel(figList)
        figH = figList(i);
        if ~isvalid(figH), continue; end

        plotKind = freq_sanitize(nameList(i));
        base = sprintf('%s__%s__%s__%02d__%s', opt.prefix, funcName, baseTitle, i, plotKind);
        base = freq_sanitize(base);
        freq_save_one_figure(figH, opt, base);
    end
end

function freq_save_one_figure(figH, opt, baseName)
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

function s = freq_numtag(x)
    % 小数点/マイナス/指数を安全にトークン化:  -0.35 -> m0p35, 1e-3 -> 1em3
    if isempty(x) || ~isfinite(x)
        s = "nan";
        return;
    end
    s = string(sprintf('%.6g', double(x)));
    s = strrep(s, "-", "m");
    s = strrep(s, "+", "");
    s = strrep(s, ".", "p");
end

function s = freq_sanitize(s)
    s = string(s);
    s = regexprep(s, '[\\/:*?"<>|]', '_');
    s = regexprep(s, '\s+', '_');
    s = regexprep(s, '_+', '_');
    s = strip(s, '_');
end
