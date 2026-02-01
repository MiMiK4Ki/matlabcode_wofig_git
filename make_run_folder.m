function outDir = make_run_folder(root, param, win, label, extra)
% extra: 任意 struct（symR, nSamps, beta, nsymb, filtergain, fIF など入れたい値）
    if nargin < 5, extra = struct(); end

    ts = datestr(now,'yyyymmdd_HHMMSS');

    modtype = "NA";
    if isfield(param,'MODTYPE') && ~isempty(param.MODTYPE)
        modtype = upper(string(param.MODTYPE));
    end

    tok = {};
    tok{end+1} = string(label);
    tok{end+1} = "mod_" + modtype;
    tok{end+1} = "M_" + string(param.MODNUM);
    tok{end+1} = "SNR_" + numtag(param.SNRdB);
    tok{end+1} = "NoSpS_" + string(param.NoSpS);
    tok{end+1} = "N_" + string(param.TXD_N);
    tok{end+1} = "bmax_" + numtag(param.bmax_initial);
    tok{end+1} = "C_" + string(win.first_C) + "to" + string(win.C_max);
    tok{end+1} = "F_" + string(win.first_F) + "to" + string(win.F_max);

    % extra も入れたいときだけ
    addNum = @(k) (isfield(extra,k) && isscalar(extra.(k)) && isfinite(extra.(k)));
    if addNum('symR'),       tok{end+1} = "symR_" + numtag(extra.symR); end
    if addNum('nSamps'),     tok{end+1} = "nSamps_" + numtag(extra.nSamps); end
    if addNum('beta'),       tok{end+1} = "beta_" + numtag(extra.beta); end
    if addNum('nsymb'),      tok{end+1} = "nsymb_" + numtag(extra.nsymb); end
    if addNum('filtergain'), tok{end+1} = "gain_" + numtag(extra.filtergain); end
    if addNum('fIF'),        tok{end+1} = "fIF_" + numtag(extra.fIF); end

    name = ts + "__" + strjoin(string(tok), "__");
    name = sanitize(name);

    outDir = fullfile(root, "results", name);
    if ~exist(outDir,'dir'), mkdir(outDir); end
end

function s = numtag(x)
    s = string(sprintf('%.6g', double(x)));
    s = strrep(s, "-", "m");
    s = strrep(s, "+", "");
    s = strrep(s, ".", "p");
end

function s = sanitize(s)
    s = string(s);
    s = regexprep(s, '[\\/:*?"<>|]', '_');
    s = regexprep(s, '\s+', '_');
    s = regexprep(s, '_+', '_');
    s = strip(s, '_');
end
