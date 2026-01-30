function gui_plot_capacity_new3()
    % GUI: Capacity Plotter V3 with clear, open window, and new model sigma plot
    % data = load(fullfile('\\NAS37092E','home','sim','SimResult34030TXDN.mat'),'param','res');
     data = load(fullfile('./SimResult20000TXDN.mat'),'param','res');
    param = data.param;
    res   = data.res;

    % Extract variables
    r_ik_cell       = res.r_ik;
    Hard_capacity   = res.Hard_capacity;
    fitResults_cell = res.fitResults;
    vardsi_list     = res.vardsi_list;   % for Var[I]

    % Parameter UI definitions
    paramInfo = {...
        struct('internal','cutoff_coeff','display','W_s'),...
        struct('internal','SNRdB','display','SNR (dB)'),...
        struct('internal','first_F','display','F'),...
        struct('internal','bmax_initial','display','b_{max}/M'),...
        struct('internal','MODNUM','display','2M')};
    nParams = numel(paramInfo);
    for i = 1:nParams
        internalNames{i} = paramInfo{i}.internal;
        displayNames{i}  = paramInfo{i}.display;
        paramValues{i}   = param.(internalNames{i});
    end

    % Create UI
    ctrlWidth = 300;
    btnColX   = 160;  utilColX = 10;   
    OFs=200;
    fig = uifigure('Name','Capacity Plotter V3','Position',[100 100 1200 700+OFs]);
    pnl = uipanel(fig,'Title','Controls','Position',[0 0 ctrlWidth 700+OFs]);
    ax  = uiaxes(fig,'Position',[ctrlWidth+20 50 850 800]); hold(ax,'on');
    ax.FontSize = 12; grid(ax,'on'); legend(ax,'off');

    % X-axis dropdown
    uilabel(pnl,'Text','X-axis:','Position',[10 640+OFs 60 22]);
    ddX = uidropdown(pnl,'Items',displayNames,'Position',[80 640+OFs 200 22],...
        'ValueChangedFcn',@(dd,~) onXChanged(dd));

    % Fixed parameter dropdowns
    y=600+OFs;
    for i=1:nParams
        uilabel(pnl,'Text',displayNames{i},'Position',[10 y 60 22]);
        dd = uidropdown(pnl,'Items',string(paramValues{i}),'Position',[80 y 200 22]);
        paramUI.(internalNames{i}) = dd;
        y=y-40;
    end

    % Plot buttons (right)
    btnY=360+OFs;
    uibutton(pnl,'Text','Plot Capacity','Position',[btnColX btnY 120 30],'ButtonPushedFcn',@(~,~) plotVar('Capacity'));
    btnY=btnY-40;
    uibutton(pnl,'Text','Plot P\_wTHP','Position',[btnColX btnY 120 30],'ButtonPushedFcn',@(~,~) plotVar('P_wTHP'));
    btnY=btnY-40;
    uibutton(pnl,'Text','Plot r\_ik','Position',[btnColX btnY 120 30],'ButtonPushedFcn',@(~,~) plotVar('r_ik'));
    btnY=btnY-40;
    uibutton(pnl,'Text','Plot mu','Position',[btnColX btnY 120 30],'ButtonPushedFcn',@(~,~) plotVar('mu'));
    btnY=btnY-40;
    uibutton(pnl,'Text','Plot sigma','Position',[btnColX btnY 120 30],'ButtonPushedFcn',@(~,~) plotVar('sigma'));
    btnY=btnY-40;
    uibutton(pnl,'Text','Plot HardCap','Position',[btnColX btnY 120 30],'ButtonPushedFcn',@(~,~) plotVar('Hard_capacity'));
    btnY=btnY-40;
    uibutton(pnl,'Text','Plot Model Sigma','Position',[btnColX btnY 120 30],'ButtonPushedFcn',@(~,~) plotVar('sigma_model'));
    btnY=btnY-40;
    uibutton(pnl,'Text','Plot Sigma Approx','Position',[btnColX btnY 120 30],'ButtonPushedFcn',@(~,~) plotVar('model_sigma_approx'));
    btnY = btnY - 40;
    uibutton(pnl, 'Text', 'Plot Var_{bk}', ...
        'Position', [btnColX, btnY, 120, 30], ...
        'ButtonPushedFcn', @(~,~) plotVar('Var_bk'));
    
    btnY = btnY - 40;
    uibutton(pnl, 'Text', 'Plot h^2_{idx>=0}', ...
        'Position', [btnColX, btnY, 120, 30], ...
        'ButtonPushedFcn', @(~,~) plotVar('sum_h2_pos'));
    
    btnY = btnY - 40;
    uibutton(pnl, 'Text', 'Plot h^2_{idx<0}', ...
        'Position', [btnColX, btnY, 120, 30], ...
        'ButtonPushedFcn', @(~,~) plotVar('sum_h2_neg'));
    
    btnY = btnY - 40;
    uibutton(pnl, 'Text', 'Plot ImpResp Power', ...
        'Position', [btnColX, btnY, 120, 30], ...
        'ButtonPushedFcn', @(~,~) plotVar('impulse_power'));    
    % Utility buttons (left)
    utilY=360;
    uibutton(pnl,'Text','Clear','Position',[utilColX utilY 120 30],'ButtonPushedFcn',@(~,~) clearPlots());
    utilY=utilY-40;
    uibutton(pnl,'Text','Open Graph Window','Position',[utilColX utilY 120 30],'ButtonPushedFcn',@(~,~) openGraphWindow());

    % Initial state
    currentX = displayNames{1}; ddX.Value = currentX;

    function onXChanged(dd)
        currentX = dd.Value; clearPlots();
    end
    function clearPlots()
        cla(ax); hold(ax,'on'); legend(ax,'off');
    end

    function plotVar(varName)
        [xVec, idxCell, fixedStr] = getXandSlice();
        iX = find(strcmp(displayNames,currentX),1);
        switch varName
            case 'Capacity'
                yv = squeeze(res.Capacity_wTHP(idxCell{:})); h=plot(ax,xVec,yv,'LineWidth',2);
            case 'P_wTHP'
                yv = squeeze(res.P_wTHP(idxCell{:})); h=plot(ax,xVec,yv,'LineWidth',2);
            case 'Hard_capacity'
                yv = squeeze(Hard_capacity(idxCell{:})); h=plot(ax,xVec,yv,'LineWidth',2);
            case 'r_ik'
                cellA = squeeze(r_ik_cell(idxCell{:})); N=numel(cellA);
                Mside = size(cellA{1},1);
                for ii=1:Mside, for kk=1:Mside
                    ytmp = arrayfun(@(j) cellA{j}(ii,kk),1:N);
                    h=plot(ax,xVec,ytmp,'LineWidth',1,'DisplayName',sprintf('r_{%d,%d}',ii,kk));
                end,end
            case 'mu'
                fr = squeeze(fitResults_cell(idxCell{:})); N=numel(fr);
                arr0 = fr{1}; S=numel(arr0);
                for s=1:S
                    sym=arr0(s).symbol;
                    ytmp = arrayfun(@(j) fr{j}(s).mu,1:N);
                    h=plot(ax,xVec,ytmp,'LineWidth',2,'DisplayName',sprintf('mu_sym%d',sym));
                end
            case 'sigma'
                fr = squeeze(fitResults_cell(idxCell{:})); N=numel(fr);
                arr0 = fr{1}; S=numel(arr0);
                for s=1:S
                    sym=arr0(s).symbol;
                    ytmp = arrayfun(@(j) fr{j}(s).sigma,1:N);
                    h=plot(ax,xVec,1./ytmp,'LineWidth',2,'DisplayName',sprintf('sigma_sym%d',sym));
                end
            case 'sigma_model'
                % Model sigma を計算するための下準備 (gui_plot_capacity_new2 を参考に)
                % P_wTHP と Normalized_Coeff
                Pvec  = squeeze(res.P_wTHP(idxCell{:}));
                NCvec = squeeze(res.Normalized_Coeff(idxCell{:}));
                % Var[I] (ikVal)
                A     = squeeze(vardsi_list(idxCell{:}));
                B     = arrayfun(@(c) c.normalized, A);
                N     = numel(B);
                % actualSNRval
                snrIdx = find(strcmp(internalNames,'SNRdB'));
                cutIdx = find(strcmp(internalNames,'cutoff_coeff'));
                actualSNRvec = zeros(1,N);
                for j = 1:N
                    if snrIdx == iX
                        snrDbVal = param.SNRdB(j);
                    else
                        snrDbVal = str2double(paramUI.SNRdB.Value);
                    end
                    snrLin = 10^(snrDbVal/10);
                    if cutIdx == iX
                        cutVal = param.cutoff_coeff(j);
                    else
                        cutVal = str2double(paramUI.cutoff_coeff.Value);
                    end
                    if cutVal >= 1
                        actualSNRvec(j) = snrLin;
                    else
                        actualSNRvec(j) = snrLin / cutVal;
                    end
                end
                % ytmp の計算
                ytmp = zeros(1,N);
                for j = 1:N
                    ikVal = B(j).ik; Pval = Pvec(j); NC = NCvec(j); snrj = actualSNRvec(j);
                    ytmp(j) = sqrt( (Pval/NC^2 )/snrj + ikVal );
                end
                h = plot(ax, xVec, ytmp, 'LineWidth', 2, 'DisplayName', 'sigma_model');
        case 'model_sigma_approx'
            % --- データ準備 ---
            % 4-D double から数値ベクトルを抜き出し
            VARBvec    = squeeze(res.VARB(idxCell{:}));             % 新しい変数 VARB
            impEvec    = squeeze(res.impulse_energ(idxCell{:}));    % impulse energy
            NCvec      = squeeze(res.Normalized_Coeff(idxCell{:})); % Normalized_Coeff
            % 4-D cell からセル配列を抜き出し
            chCcell    = squeeze(res.ChannelCoeff(idxCell{:}));     % ChannelCoeff cell array

            % インデックス取得
            snrIdx    = find(strcmp(internalNames,'SNRdB'));
            cutIdx    = find(strcmp(internalNames,'cutoff_coeff'));
            firstFIdx = find(strcmp(internalNames,'first_F'));
            modIdx    = find(strcmp(internalNames,'MODNUM'));

            N = numel(VARBvec);
            actualSNRvec = zeros(1,N);
            cutoffVec    = zeros(1,N);
            firstFvec    = zeros(1,N);
            modNumVec    = zeros(1,N);
            for j = 1:N
                % actual SNR
                if snrIdx == iX
                    snrDbVal = param.SNRdB(j);
                else
                    snrDbVal = str2double(paramUI.SNRdB.Value);
                end
                snrLin = 10^(snrDbVal/10);

                % cutoff_coeff
                if cutIdx == iX
                    cutVal = param.cutoff_coeff(j);
                else
                    cutVal = str2double(paramUI.cutoff_coeff.Value);
                end
                cutoffVec(j) = cutVal;
                if cutVal >= 1
                    actualSNRvec(j) = snrLin;
                else
                    actualSNRvec(j) = snrLin / cutVal;
                end

                % first_F
                if firstFIdx == iX
                    firstFvec(j) = param.first_F(j);
                else
                    firstFvec(j) = str2double(paramUI.first_F.Value);
                end

                % MODNUM
                if modIdx == iX
                    modNumVec(j) = param.MODNUM(j);
                else
                    modNumVec(j) = str2double(paramUI.MODNUM.Value);
                end
            end

            % --- 近似モデル σ の計算用ループ ---
            ytmp = zeros(1,N);
            for j = 1:N
                varb      = VARBvec(j);               % 新しい変数 VARB
                impE      = impEvec(j);               % impulse energy
                NC        = NCvec(j);                 % Normalized_Coeff
                chC       = chCcell{j};               % Channel_Coeff_scalar
                snrj      = actualSNRvec(j);          % actual SNR
                Tc_val    = param.Tc_but / cutoffVec(j);
                firstFval = firstFvec(j);
                Mval      = modNumVec(j);

                % ここに近似モデルの数式を記述してください:
                first_C = -2;
                VARB=Mval^2/3
                % VARB=(Mval^2-1)/3 / sum( chC((1+firstFval-first_C):end).^2)
                VARB = varb
                ytmp(j) = sqrt(VARB * ( 1/snrj/Tc_val *impE/(NC^2) + sum( chC(1:(firstFval-first_C)).^2) ) ) ;
            end

            % プロット
            h = plot(ax, xVec, ytmp, 'LineWidth', 2, 'DisplayName','sigma_model_approx');
        case 'Var_bk'
            % res.VARB は 4-D double → ベクトル化してそのままプロット
            yv = squeeze(res.VARB(idxCell{:}));
            h = plot(ax, xVec, yv, 'LineWidth',2, 'DisplayName','Var_{bk}');

        case 'sum_h2_pos'
            % h^2 の 0 以上インデックス和
            chCcell = squeeze(res.ChannelCoeff(idxCell{:}));  % 1×N cell
            snrIdx    = find(strcmp(internalNames,'SNRdB'));
            cutIdx    = find(strcmp(internalNames,'cutoff_coeff'));
            firstFIdx = find(strcmp(internalNames,'first_F'));
            N = numel(chCcell);
            ytmp = zeros(1,N);
            first_C = -2;  % システム定数
            for j=1:N
                % first_F の取得
                if firstFIdx==iX
                    firstFval = param.first_F(j);
                else
                    firstFval = str2double(paramUI.first_F.Value);
                end
                % チャネル係数 h
                hvec = chCcell{j};
                h2   = hvec.^2;
                startIdx = 1 + firstFval - first_C;
                ytmp(j) = sum( h2(startIdx:end) );
            end
            h = plot(ax, xVec, ytmp, 'LineWidth',2, 'DisplayName','sum\,h^2_{idx\ge0}');

        case 'sum_h2_neg'
            % h^2 の負インデックス和
            chCcell = squeeze(res.ChannelCoeff(idxCell{:}));
            firstFIdx = find(strcmp(internalNames,'first_F'));
            N = numel(chCcell);
            ytmp = zeros(1,N);
            first_C = -2;
            for j=1:N
                if firstFIdx==iX
                    firstFval = param.first_F(j);
                else
                    firstFval = str2double(paramUI.first_F.Value);
                end
                hvec = chCcell{j};
                h2   = hvec.^2;
                endIdx = firstFval - first_C;
                ytmp(j) = sum( h2(1:endIdx) );
            end
            h = plot(ax, xVec, ytmp, 'LineWidth',2, 'DisplayName','sum\,h^2_{idx<0}');

        case 'impulse_power'
            % 1/Tc_val * impE/(NC^2)
            impEvec = squeeze(res.impulse_energ(idxCell{:}));    % 4-D double
            NCvec   = squeeze(res.Normalized_Coeff(idxCell{:})); % 4-D double
            cutIdx  = find(strcmp(internalNames,'cutoff_coeff'));
            N = numel(impEvec);
            ytmp = zeros(1,N);
            for j=1:N
                % cutoff_coeff の取得
                if cutIdx==iX
                    cutVal = param.cutoff_coeff(j);
                else
                    cutVal = str2double(paramUI.cutoff_coeff.Value);
                end
                Tc_val = param.Tc_but / cutVal;
                ytmp(j) = (1/Tc_val) * ( impEvec(j) / (NCvec(j)^2) );
            end
            h = plot(ax, xVec, ytmp, 'LineWidth',2, 'DisplayName','ImpResp\,Power');

                
            otherwise
                return;
        end
        if ~isempty(fixedStr)
            h.DisplayName = [h.DisplayName ' (' fixedStr ')'];
        end
        legend(ax,'show','Location','best');
    end

    function [xVec, idxCell, fixedStr] = getXandSlice()
        iX = find(strcmp(displayNames,currentX),1);
        xVec = param.(internalNames{iX});
        idxCell = cell(1,nParams);
        fixedStr = '';
        for ii=1:nParams
            if ii==iX, idxCell{ii}=1:numel(paramValues{ii});
            else
                val = str2double(paramUI.(internalNames{ii}).Value);
                idx = find(abs(paramValues{ii}-val)<1e-8,1);
                idxCell{ii}=idx;
                fixedStr = [fixedStr displayNames{ii} '=' num2str(val) ', '];
            end
        end
        if ~isempty(fixedStr), fixedStr(end-1:end)=[]; end
    end

    function openGraphWindow()
        fig2=figure('Name','Graph Copy','NumberTitle','off');
        ax2=axes('Parent',fig2);
        copyobj(allchild(ax),ax2);
        ax2.Title.String=ax.Title.String;
        ax2.XLabel.String=ax.XLabel.String;
        ax2.YLabel.String=ax.YLabel.String;
        legend(ax2,'show','Location','best'); grid(ax2,'on');
    end
end
