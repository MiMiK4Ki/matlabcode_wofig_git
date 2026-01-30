function gui_plot_sweep3routes(matFile)
% gui_plot_sweep3routes
%   run_sweep_3routes_parfor.m が保存した
%     - paramSweep
%     - res
%   を読み込み、X軸=走査パラメータ、Y軸=任意の結果をGUIでプロットする。
%
% 追加:
%   - Max SoftCap * W_s
%   - Max HardCap * W_s
%   → X軸の各点ごとに、他パラメータ全探索での最大値を描く
%
% 使い方:
%   gui_plot_sweep3routes();
%   gui_plot_sweep3routes('Sweep3routes_20000sym.mat');

if nargin < 1 || isempty(matFile)
    matFile = 'Sweep3routes_20050sym.mat';
    if ~exist(matFile,'file')
        [f,p] = uigetfile('*.mat','Select sweep result .mat');
        if isequal(f,0), return; end
        matFile = fullfile(p,f);
    end
end

data = load(matFile);

if isfield(data,'paramSweep')
    param = data.paramSweep;
elseif isfield(data,'param')
    param = data.param;
else
    error('This mat file has no paramSweep/param.');
end
if isfield(data,'res')
    res = data.res;
else
    error('This mat file has no res.');
end

% -------------------------
% Sweep次元（gridboxの順序が前提）
% gridbox = [cutoff, SNR, first_F, bmax, MODNUM]
% -------------------------
paramInfo = { ...
    struct('internal','cutoff_coeff','display','W_s (=cutoff)'), ...
    struct('internal','SNRdB','display','SNR (dB)'), ...
    struct('internal','first_F','display','F'), ...
    struct('internal','bmax_initial','display','b_{max}/M'), ...
    struct('internal','MODNUM','display','2M') };

nParams = numel(paramInfo);
internalNames = cell(1,nParams);
displayNames  = cell(1,nParams);
paramValues   = cell(1,nParams);
for i = 1:nParams
    internalNames{i} = paramInfo{i}.internal;
    displayNames{i}  = paramInfo{i}.display;
    paramValues{i}   = param.(internalNames{i});
end

% gridサイズ
nC = numel(param.cutoff_coeff);
nS = numel(param.SNRdB);
nF = numel(param.first_F);
nB = numel(param.bmax_initial);
nM = numel(param.MODNUM);

% cutoff の5次元グリッド（SoftCap*W_s 等を「正しく」作るため）
% サイズ = [nC nS nF nB nM]
CUTGRID = repmat(reshape(param.cutoff_coeff(:), [nC 1 1 1 1]), [1 nS nF nB nM]);

% -------------------------
% GUI
% -------------------------
ctrlWidth = 360;
fig = uifigure('Name',['Sweep Plotter: ' matFile],'Position',[100 100 1280 840]);
pnl = uipanel(fig,'Title','Controls','Position',[0 0 ctrlWidth 840]);
ax  = uiaxes(fig,'Position',[ctrlWidth+20 60 880 750]);
hold(ax,'on'); grid(ax,'on'); legend(ax,'off');
ax.FontSize = 12;

% X-axis dropdown
uilabel(pnl,'Text','X-axis:','Position',[10 790 60 22]);
ddX = uidropdown(pnl,'Items',displayNames,'Position',[80 790 260 22]);

% 固定パラメータ dropdown
paramUI = struct();
y0 = 750;
for i = 1:nParams
    uilabel(pnl,'Text',displayNames{i},'Position',[10 y0 110 22]);
    dd = uidropdown(pnl,'Items',string(paramValues{i}),'Position',[130 y0 210 22]);
    paramUI.(internalNames{i}) = dd;
    y0 = y0 - 40;
end

% Route dropdown
uilabel(pnl,'Text','Route:','Position',[10 550 60 22]);
% Route dropdown
uilabel(pnl,'Text','Route:','Position',[10 550 60 22]);
ddRoute = uidropdown(pnl,'Items',{ ...
    'Theory (th)', ...
    'PartialExp (px)', ...
    'S2P (s2p)', ...
    'Both (th+px)', ...
    'Both (th+s2p)', ...
    'Both (px+s2p)', ...
    'All (th+px+s2p)'}, ...
    'Position',[80 550 260 22], 'Value','Theory (th)');


% mode dropdown
uilabel(pnl,'Text','Mode:','Position',[10 510 60 22]);
ddMode = uidropdown(pnl,'Items',{'woTHP','wTHP'}, ...
    'Position',[80 510 260 22], 'Value','woTHP');

% Y variable dropdown
uilabel(pnl,'Text','Y-axis:','Position',[10 470 60 22]);
ddY = uidropdown(pnl,'Items',{ ...
    'BER', ...
    'HardCap', 'HardCap * W_s', ...
    'SoftCap', 'SoftCap * W_s', ...
    'EVM', 'EVMref', ...
    'Max SoftCap * W_s', ...
    'Max HardCap * W_s' ...
    }, ...
    'Position',[80 470 260 22], 'Value','SoftCap');


% Buttons
uibutton(pnl,'Text','Plot','Position',[220 420 120 32], ...
    'ButtonPushedFcn',@(~,~) plotSelected());
uibutton(pnl,'Text','Clear','Position',[90 420 120 32], ...
    'ButtonPushedFcn',@(~,~) clearPlots());
uibutton(pnl,'Text','Open Graph Window','Position',[90 375 250 32], ...
    'ButtonPushedFcn',@(~,~) openGraphWindow());

ddX.ValueChangedFcn = @(~,~) clearPlots();

% -------------------------
% callbacks
% -------------------------
    function clearPlots()
        cla(ax);
        hold(ax,'on'); grid(ax,'on');
        legend(ax,'off');
        set(ax,'YScale','linear');
    end

    function plotSelected()
        [xVec, idxCellFixed, fixedStr, iX] = getXandSlice();
        routeSel = ddRoute.Value;
        modeSel  = ddMode.Value;
        ySel     = ddY.Value;

        % route選択
        % route選択
        routeKeys = {};
        routeLabels = {};

        switch ddRoute.Value
            case 'Theory (th)'
                routeKeys = {'th'};     routeLabels = {'th'};
            case 'PartialExp (px)'
                routeKeys = {'px'};     routeLabels = {'px'};
            case 'S2P (s2p)'
                routeKeys = {'s2p'};    routeLabels = {'s2p'};
            case 'Both (th+px)'
                routeKeys = {'th','px'};      routeLabels = {'th','px'};
            case 'Both (th+s2p)'
                routeKeys = {'th','s2p'};     routeLabels = {'th','s2p'};
            case 'Both (px+s2p)'
                routeKeys = {'px','s2p'};     routeLabels = {'px','s2p'};
            case 'All (th+px+s2p)'
                routeKeys = {'th','px','s2p'}; routeLabels = {'th','px','s2p'};
            otherwise
                error('Unknown route selection: %s', ddRoute.Value);
        end


        for ir = 1:numel(routeKeys)
            routeKey = routeKeys{ir};

            [baseField, yLabel, needCutoffScale, isMaxMode] = mapField(routeKey, ySel, modeSel);

            if isMaxMode
                % ======================================================
                % Maxモード:
                %   X軸の各点ごとに、他パラメータ全探索で max を取る
                %   (SoftCap*W_s 等は「掛けてから」max)
                % ======================================================
                yv = computeMaxCurve(baseField, iX, needCutoffScale);

                h = plot(ax, xVec, yv, 'LineWidth',2);
                h.DisplayName = sprintf('%s %s %s (MAX over others)', routeLabels{ir}, ySel, modeSel);

                xlabel(ax, ddX.Value);
                ylabel(ax, yLabel);
                set(ax,'YScale','linear');  % Max系は基本linear

            else
                % ======================================================
                % 通常モード（固定条件で断面をプロット）
                % ======================================================
                yv = squeeze(res.(baseField)(idxCellFixed{:}));
                yv = yv(:);

                if needCutoffScale
                    % cutoff を「固定条件に応じて」掛ける
                    cutIdx = find(strcmp(internalNames,'cutoff_coeff'),1);
                    if iX == cutIdx
                        cutVec = xVec(:);
                    else
                        cutFixed = str2double(paramUI.cutoff_coeff.Value);
                        cutVec = cutFixed * ones(numel(xVec),1);
                    end
                    yv = yv .* cutVec;
                end

                h = plot(ax, xVec, yv, 'LineWidth',2);
                h.DisplayName = sprintf('%s %s %s (%s)', routeLabels{ir}, ySel, modeSel, fixedStr);

                xlabel(ax, ddX.Value);
                ylabel(ax, yLabel);

                if strcmp(ySel,'BER')
                    set(ax,'YScale','log');
                else
                    set(ax,'YScale','linear');
                end
            end
        end

        legend(ax,'show','Location','best');
        title(ax, sprintf('%s | X=%s', matFile, ddX.Value));
    end

    function yv = computeMaxCurve(fieldName, iX, needCutoffScale)
        % X軸の各点ごとに、他次元を全探索して最大値を返す
        nX = numel(paramValues{iX});
        yv = NaN(nX,1);

        for ix = 1:nX
            idxCell = cell(1,nParams);
            for d = 1:nParams
                if d == iX
                    idxCell{d} = ix;
                else
                    idxCell{d} = 1:numel(paramValues{d});  % 全探索
                end
            end

            sub = res.(fieldName)(idxCell{:});

            if needCutoffScale
                sub = sub .* CUTGRID(idxCell{:});  % (metric * cutoff) を作ってから max
            end

            v = sub(:);
            v = v(~isnan(v));   % pxがNaN混じりでも最大が取れるように
            if ~isempty(v)
                yv(ix) = max(v);
            end
        end
    end

    function [fieldName, yLabel, needCutoffScale, isMaxMode] = mapField(routeKey, ySel, modeSel)
        needCutoffScale = false;
        isMaxMode = false;

        switch ySel
            case 'BER'
                metric = 'BER';
                yLabel = 'BER';
            case 'HardCap'
                metric = 'HardCap';
                yLabel = 'Hard Capacity';
            case 'HardCap * W_s'
                metric = 'HardCap';
                yLabel = 'Hard Capacity * W_s';
                needCutoffScale = true;

            case 'SoftCap'
                metric = 'SoftCap';
                yLabel = 'Soft Capacity';
            case 'SoftCap * W_s'
                metric = 'SoftCap';
                yLabel = 'Soft Capacity * W_s';
                needCutoffScale = true;

            case 'EVM'
                metric = 'EVM';
                yLabel = 'EVM (%)';

            case 'Max SoftCap * W_s'
                metric = 'SoftCap';
                yLabel = 'max(SoftCap * W_s)';
                needCutoffScale = true;
                isMaxMode = true;

            case 'Max HardCap * W_s'
                metric = 'HardCap';
                yLabel = 'max(HardCap * W_s)';
                needCutoffScale = true;
                isMaxMode = true;
            case 'EVMref'
                metric = 'EVMref';
                yLabel = 'EVMref (%)';

            otherwise
                error('Unknown Y selection: %s', ySel);
        end

        fieldName = sprintf('%s_%s_%s', routeKey, metric, modeSel);

        if ~isfield(res, fieldName)
            error('Field not found in res: %s', fieldName);
        end
    end

    function [xVec, idxCell, fixedStr, iX] = getXandSlice()
        iX = find(strcmp(displayNames, ddX.Value), 1);
        xVec = param.(internalNames{iX});
        xVec = xVec(:);

        idxCell = cell(1, nParams);
        fixedStr = '';

        for ii = 1:nParams
            if ii == iX
                idxCell{ii} = 1:numel(paramValues{ii});
            else
                val = str2double(paramUI.(internalNames{ii}).Value);
                vlist = paramValues{ii}(:).';
                idx = find(abs(vlist - val) < 1e-12, 1);
                idxCell{ii} = idx;
                fixedStr = [fixedStr, displayNames{ii}, '=', num2str(val), ', ']; %#ok<AGROW>
            end
        end
        if ~isempty(fixedStr)
            fixedStr(end-1:end) = [];
        end
    end

    function openGraphWindow()
        fig2 = figure('Name','Graph Copy','NumberTitle','off');
        ax2  = axes('Parent',fig2); hold(ax2,'on'); grid(ax2,'on');

        ch = ax.Children;
        for k = numel(ch):-1:1
            try
                copyobj(ch(k), ax2);
            catch
            end
        end
        xlabel(ax2, ax.XLabel.String);
        ylabel(ax2, ax.YLabel.String);
        title(ax2, ax.Title.String);
        legend(ax2,'show','Location','best');
        ax2.FontSize=28;
    end
end
