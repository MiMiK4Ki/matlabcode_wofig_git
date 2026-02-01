%% ================= main_loadplot_saved_results_by_date.m =================
% 目的:
%   - yyyymmdd で保存された "runフォルダ" 群から workspace .mat を読み込み、
%     指定した横軸(X)に対して、指定した縦軸(Y)指標(EVM/Capacity/BER/SER 等)をプロットする。
%
% 前提(保存側):
%   - runフォルダ名が "yyyymmdd__..." で始まる（ts=datestr(now,'yyyymmdd') を使用）
%   - 各 runフォルダ直下に workspace mat (例: workspace_all.mat) が保存されている
%
% 使い方:
%   1) 下の「USER SETTINGS」を編集
%   2) このスクリプトを実行
%
% 追加したい軸の候補:
%   - X軸: get_x_value() に case を追加
%   - Y軸: get_metric_value() に case を追加
%
% NOTE:
%   - 同じX値が複数runに存在する場合は doGroupAverage=true で平均化してプロット可能

clear; close all; clc;

%% ================= USER SETTINGS =================
% (1) 読み込み対象の日付 (yyyymmdd) を手入力
dateTag = "20260201";   % ★ここを毎回変更

% (2) runフォルダの親ディレクトリ
% 例: 保存側で outDir = fullfile(pwd,'results', name) にしているなら "results"
rootDir = fullfile(pwd, "results");   % ★保存側に合わせて変更

% (3) 各 runフォルダ内の workspace mat ファイル名
matFileName = "workspace_all.mat";    % ★保存側の名前に合わせる

% (4) フォルダ名の部分一致フィルタ（不要なら ""）
% 例: "Dsimple_v2" 等を入れると、その文字列を含むrunだけ拾う
folderMustContain = "";

% (5) 横軸（切り替えたいものを1つ指定）
% 例: "MODNUM", "symR", "F", "bmax", "SNRdB"
xKey = "SNRdB";

% (6) 縦軸（複数指定OK）
% "capacity" は Res.Soft_capacity を指す（あなたの evaluate 系の命名に合わせる）
yKeys = ["EVM","capacity","hard_capacity","BER","SER"];

% (7) 表示するモード（保存matに入っているものだけ出る）
% 典型: Dルートなら woTHP=Res_woE, wTHP=Res_wE
modes = ["woTHP","wTHP"];

% (8) 同一X値が複数runあるときの平均化
doGroupAverage = true;

% (9) BER/SER をログ軸にするか
useLogY_for_BER_SER = true;

% (10) 固定条件フィルタ（不要なら空のstructでOK）
% mainの param 定義のように並べて「一致したrunだけ」を拾う用途
%   - 文字列: "" (空) は無視
%   - 数値: NaN は無視
fixed = struct();
fixed.MODTYPE = "QAM";     % "PAM" / "QAM" / ""(無視)
fixed.SNRdB   = NaN;    % 例: 5
fixed.NoSpS   = 50;    % 例: 5
fixed.first_F = 0;    % 例: -1
fixed.symR    = 100e6;    % 例: 100e6
fixed.bmax_initial = 2;
fixed.MODNUM  = 4;

%% ================= load all runs in the day =================
runDirs = dir(char(fullfile(rootDir, dateTag + "__*")));
runDirs = runDirs([runDirs.isdir]);

if ~isempty(folderMustContain)
    keep = false(size(runDirs));
    for i = 1:numel(runDirs)
        keep(i) = contains(string(runDirs(i).name), string(folderMustContain));
    end
    runDirs = runDirs(keep);
end

if isempty(runDirs)
    error('No run folders found. rootDir=%s, dateTag=%s', rootDir, dateTag);
end

fprintf('Found %d run folders under %s for date %s\n', numel(runDirs), rootDir, dateTag);

runs = struct([]);
for i = 1:numel(runDirs)
    folderPath = fullfile(runDirs(i).folder, runDirs(i).name);

    try
        R = load_one_run(folderPath, matFileName);
    catch ME
        warning('[skip] %s : load failed (%s)', runDirs(i).name, ME.message);
        continue;
    end

    % fixed filter
    if ~match_fixed(R, fixed)
        continue;
    end

    runs = [runs; R]; %#ok<AGROW>
end

if isempty(runs)
    error('No runs left after filtering. Please relax fixed/folderMustContain.');
end

fprintf('Loaded %d runs after filtering.\n', numel(runs));

%% ================= build x / y arrays =================
xVals = nan(numel(runs),1);
for i = 1:numel(runs)
    xVals(i) = get_x_value(runs(i), xKey);
end

% drop runs where x is NaN
okX = isfinite(xVals);
runs = runs(okX);
xVals = xVals(okX);

if isempty(runs)
    error('All runs dropped because X="%s" could not be extracted.', xKey);
end

% x scaling/label
[xPlot, xLabel] = x_scale_and_label(xVals, xKey);

%% ================= plot per metric =================
for yk = 1:numel(yKeys)
    yKey = string(yKeys(yk));

    % collect y per mode
    Y = nan(numel(runs), numel(modes));
    for mi = 1:numel(modes)
        md = string(modes(mi));
        for i = 1:numel(runs)
            Y(i,mi) = get_metric_value(runs(i), md, yKey);
        end
    end

    % sorting by x
    [xPlot_s, idx] = sort(xPlot);
    xVals_s = xVals(idx);
    Y_s = Y(idx,:);

    if doGroupAverage
        [xU, YU, Ystd] = group_average(xVals_s, xPlot_s, Y_s);
        x_use = xU;
        Y_use = YU;
        Yerr  = Ystd;
    else
        x_use = xPlot_s;
        Y_use = Y_s;
        Yerr  = [];
    end

    % plot
    fig = figure; hold on; grid on;

    % logY for BER/SER
    doLogY = false;
    if useLogY_for_BER_SER
        if any(strcmpi(yKey, ["BER","SER"]))
            doLogY = true;
        end
    end

    markers = ["o","s","^","d","x","+","*","v",">","<"];
    for mi = 1:numel(modes)
        mk = markers(1+mod(mi-1, numel(markers)));
        yv = Y_use(:,mi);

        if doLogY
            semilogy(x_use, yv, ['-' char(mk)], 'LineWidth', 1.5, 'DisplayName', char(modes(mi)));
        else
            plot(x_use, yv, ['-' char(mk)], 'LineWidth', 1.5, 'DisplayName', char(modes(mi)));
        end

        % errorbar (平均化している場合のみ)
        if doGroupAverage && ~isempty(Yerr)
            ysd = Yerr(:,mi);
            % BER/SER のログ軸では errorbar が見づらいので、基本は描かない（必要ならここを編集）
            if ~doLogY
                try
                    errorbar(x_use, yv, ysd, 'LineStyle','none', 'HandleVisibility','off');
                catch
                end
            end
        end
    end

    xlabel(xLabel);
    ylabel(y_label(yKey));
    title(sprintf('%s vs %s  (date=%s, Nrun=%d)', y_label(yKey), xKey, dateTag, numel(runs)));
    legend('show','Location','best');

    % Xが整数のときは見やすく
    if any(strcmpi(xKey, ["MODNUM","F","first_F"]))
        xticks(unique(x_use));
    end
end

%% ================= optional: show summary table in workspace =================
summary = make_summary_table(runs, xKey, xVals, modes, yKeys);
disp(summary);

%% ================= end =================



%% ========================================================================
%% local functions
%% ========================================================================

function R = load_one_run(folderPath, matFileName)
% 1 runフォルダから workspace mat を読み込み、必要情報を抽出して R を返す

matPath = fullfile(folderPath, matFileName);
if ~exist(matPath,'file')
    % fallback: 探す（workspace/ all を含む mat を優先）
    mats = dir(fullfile(folderPath, "*.mat"));
    if isempty(mats)
        error('No .mat found in %s', folderPath);
    end
    names = string({mats.name});
    idx = find(contains(lower(names), "workspace"), 1);
    if isempty(idx), idx = find(contains(lower(names), "all"), 1); end
    if isempty(idx), idx = 1; end
    matPath = fullfile(mats(idx).folder, mats(idx).name);
end

S = load(matPath);

% param を拾う（候補の優先順）
param = pick_first_struct(S, ["paramD","param_exp","param_s2p","param"]);
if isempty(param)
    param = struct();
end

% win もあれば拾う
win = pick_first_struct(S, ["win","winD"]);
if isempty(win)
    win = struct();
end

% extra: symR, nSamps など
extra = struct();
extra.symR   = pick_first_scalar(S, ["symR","symR_exp"]);
extra.nSamps = pick_first_scalar(S, ["nSamps","NoSpS","NoSpS_exp"]);
extra.fIF    = pick_first_scalar(S, ["fIF","fc","param_fc"]);
extra.folder = string(folderPath);
extra.mat    = string(matPath);

% results: modeごとに拾う（候補の優先順は必要に応じて追加）
Res = struct();
Res.woTHP = pick_first_struct(S, ["Res_woE","Res_th_wo","Res_px_wo","Res_s2p_wo","Res_wo"]);
Res.wTHP  = pick_first_struct(S, ["Res_wE","Res_th_w","Res_px_w","Res_s2p_w","Res_w"]);

R = struct();
R.folderPath = string(folderPath);
R.folderName = string(get_last_folder(folderPath));
R.matPath    = string(matPath);
R.param      = param;
R.win        = win;
R.extra      = extra;
R.Res        = Res;
end

function name = get_last_folder(p)
    p = char(p);
    [~,name] = fileparts(p);
end

function st = pick_first_struct(S, names)
    st = [];
    for i = 1:numel(names)
        nm = char(names(i));
        if isfield(S, nm) && isstruct(S.(nm))
            st = S.(nm);
            return;
        end
    end
end

function v = pick_first_scalar(S, names)
    v = NaN;
    for i = 1:numel(names)
        nm = char(names(i));
        if isfield(S, nm)
            tmp = S.(nm);
            if isnumeric(tmp) && isscalar(tmp)
                v = double(tmp);
                return;
            end
        end
    end
end

function tf = match_fixed(R, fixed)
% fixed に指定した条件を満たすrunだけ通す
% ルール:
%   - fieldが無い / 空文字 / NaN は無視
    tf = true;
    if isempty(fixed) || ~isstruct(fixed)
        return;
    end

    fns = fieldnames(fixed);
    for i = 1:numel(fns)
        key = fns{i};
        val = fixed.(key);

        % ignore rules
        if isstring(val) || ischar(val)
            if strlength(string(val)) == 0
                continue;
            end
        elseif isnumeric(val)
            if isempty(val) || (isscalar(val) && isnan(val))
                continue;
            end
        else
            continue;
        end

        rv = get_value_for_key(R, key);
        if ~isfinite_numeric_or_string(rv)
            tf = false; return;
        end

        if isstring(val) || ischar(val)
            tf = strcmpi(string(rv), string(val));
        else
            % numeric compare with tolerance
            tf = numeric_equal(double(rv), double(val));
        end

        if ~tf
            return;
        end
    end
end

function ok = isfinite_numeric_or_string(v)
    if isstring(v) || ischar(v)
        ok = (strlength(string(v)) > 0);
    elseif isnumeric(v)
        ok = isfinite(double(v));
    else
        ok = false;
    end
end

function ok = numeric_equal(a,b)
    if ~isfinite(a) || ~isfinite(b), ok = false; return; end
    if abs(round(a)-a) < 1e-12 && abs(round(b)-b) < 1e-12
        ok = (a == b);
        return;
    end
    tol = 1e-9;
    ok = (abs(a-b) <= tol*(1+abs(b)));
end

function rv = get_value_for_key(R, key)
% fixed用の値取得（param/extra/win から探す）
    key = char(key);

    if isfield(R.param, key), rv = R.param.(key); return; end
    if isfield(R.extra, key), rv = R.extra.(key); return; end
    if isfield(R.win,   key), rv = R.win.(key);   return; end

    % 別名サポート
    switch lower(key)
        case 'modnum'
            if isfield(R.param,'MODNUM'), rv = R.param.MODNUM; else, rv = NaN; end
        case 'modtype'
            if isfield(R.param,'MODTYPE'), rv = string(R.param.MODTYPE); else, rv = ""; end
        case 'snrdb'
            if isfield(R.param,'SNRdB'), rv = R.param.SNRdB; else, rv = NaN; end
        case 'first_f'
            if isfield(R.param,'first_F'), rv = R.param.first_F; else, rv = NaN; end
        case 'symr'
            rv = R.extra.symR;
        otherwise
            rv = NaN;
    end
end

function x = get_x_value(R, xKey)
% X軸にしたい値を run から取り出す
% 追加したい場合は case を増やす
    xKey = lower(string(xKey));
    p = R.param;

    switch xKey
        case "modnum"
            if isfield(p,'MODNUM'), x = double(p.MODNUM); else, x = NaN; end

        case "symr"
            x = double(R.extra.symR);

        case {"f","first_f"}
            if isfield(p,'first_F')
                x = double(p.first_F);
            elseif isfield(R.win,'first_F')
                x = double(R.win.first_F);
            else
                x = NaN;
            end

        case "bmax"
            x = calc_bmax_from_param(p);

        case "bmax_initial"
            if isfield(p,'bmax_initial'), x = double(p.bmax_initial); else, x = NaN; end
        case "SNR"
            if isfield(p,'SNRdB'), x = double(p.SNRdB); else, x = NaN; end

        otherwise
            % paramのフィールド名と同じならそれを拾う
            if isfield(p, char(xKey))
                tmp = p.(char(xKey));
                if isnumeric(tmp) && isscalar(tmp)
                    x = double(tmp);
                else
                    x = NaN;
                end
            else
                x = NaN;
            end
    end
end

function b = calc_bmax_from_param(p)
% thp_prepare_tx と同じ定義で "実際の bmax" を返す
%   - PAM: bmax = bmax_initial * (M/2)
%   - QAM: bmax = bmax_initial * (sqrt(M)/2)
    if ~isfield(p,'MODNUM') || ~isfield(p,'bmax_initial')
        b = NaN; return;
    end
    M = double(p.MODNUM);
    b0 = double(p.bmax_initial);

    modtype = "PAM";
    if isfield(p,'MODTYPE') && ~isempty(p.MODTYPE)
        modtype = upper(string(p.MODTYPE));
    end

    if modtype == "QAM"
        L = sqrt(M);
        if abs(L - round(L)) > 1e-12
            b = NaN; return;
        end
        b = b0 * (L/2);
    else
        b = b0 * (M/2);
    end
end

function y = get_metric_value(R, mode, yKey)
% Y軸(指標)を run から取り出す
% 追加したい場合は get_metric_value の case を増やす

    mode = string(mode);
    yKey = lower(string(yKey));

    if ~isfield(R.Res, char(mode)) || isempty(R.Res.(char(mode)))
        y = NaN; return;
    end
    Res = R.Res.(char(mode));

    y = metric_from_res(Res, yKey);
end

function y = metric_from_res(Res, yKey)
% Res 構造体から指標を取る
% 追加したいときはここに case を増やす
    switch yKey
        case "evm"
            y = getfield_safe(Res, 'EVM');

        case "evmref"
            y = getfield_safe(Res, 'EVMref');

        case {"capacity","soft_capacity","softcap"}
            y = getfield_safe(Res, 'Soft_capacity');

        case {"hard_capacity","hardcap"}
            y = getfield_safe(Res, 'Hard_capacity');

        case "ber"
            y = getfield_safe(Res, 'BER');

        case "ser"
            y = getfield_safe(Res, 'SER');

        otherwise
            % Resのフィールド名を直指定したい場合（例: "shift_sym"）
            fn = char(yKey);
            y = getfield_safe(Res, fn);
    end
end

function v = getfield_safe(S, fn)
    v = NaN;
    if isstruct(S) && isfield(S, fn)
        tmp = S.(fn);
        if isnumeric(tmp) && isscalar(tmp)
            v = double(tmp);
        end
    end
end

function [xPlot, xLabel] = x_scale_and_label(xVals, xKey)
% 見やすさのためのスケール変換（必要ならここを拡張）
    xKey = lower(string(xKey));
    switch xKey
        case "symr"
            xPlot = xVals / 1e6;
            xLabel = 'symR [MHz]';
        otherwise
            xPlot = xVals;
            xLabel = char(xKey);
    end
end

function yl = y_label(yKey)
% 表示用ラベル（必要なら拡張）
    yKey = lower(string(yKey));
    switch yKey
        case "evm"
            yl = 'EVM [%]';
        case "evmref"
            yl = 'EVMref [%]';
        case {"capacity","soft_capacity","softcap"}
            yl = 'Capacity (Soft) [bits/use]';
        case {"hard_capacity","hardcap"}
            yl = 'Capacity (Hard) [bits/use]';
        case "ber"
            yl = 'BER';
        case "ser"
            yl = 'SER';
        otherwise
            yl = char(yKey);
    end
end

function [xU_plot, YU, Ystd] = group_average(xRaw, xPlot, Y)
% xRaw: 元の物理量 (例: symR[Hz]) で grouping
% xPlot: 表示用スケール値
    [xU_raw, ~, gi] = unique(xRaw, 'stable');
    xU_plot = zeros(numel(xU_raw),1);
    YU   = nan(numel(xU_raw), size(Y,2));
    Ystd = nan(numel(xU_raw), size(Y,2));

    for k = 1:numel(xU_raw)
        idx = (gi == k);
        xU_plot(k) = mean(xPlot(idx));

        for c = 1:size(Y,2)
            yv = Y(idx,c);
            yv = yv(isfinite(yv));
            if isempty(yv)
                YU(k,c) = NaN;
                Ystd(k,c) = NaN;
            else
                YU(k,c) = mean(yv);
                if numel(yv) >= 2
                    Ystd(k,c) = std(yv);
                else
                    Ystd(k,c) = 0;
                end
            end
        end
    end
end

function T = make_summary_table(runs, xKey, xVals, modes, yKeys)
% 結果を table にまとめて、後で再利用しやすくする
    N = numel(runs);
    folder = strings(N,1);
    mat    = strings(N,1);
    xcol   = xVals(:);

    for i = 1:N
        folder(i) = runs(i).folderName;
        mat(i) = runs(i).matPath;
    end

    T = table(folder, mat, xcol, 'VariableNames', {'folder','mat','x'});
    T.Properties.Description = sprintf('summary: xKey=%s', xKey);

    for mi = 1:numel(modes)
        md = string(modes(mi));
        for yi = 1:numel(yKeys)
            yk = string(yKeys(yi));
            col = nan(N,1);
            for i = 1:N
                col(i) = get_metric_value(runs(i), md, yk);
            end
            vn = matlab.lang.makeValidName(md + "_" + yk);
            T.(vn) = col;
        end
    end
end
