function fitResults = fitFoldedGaussianPerSymbol2(ck_wnoise, PAM_signal, first_C, first_F, M,Mk)
% fitFoldedGaussianPerSymbol
%  シンボルごとに ck_wnoise をトランケートし、
%  生データにガウスフィッティングおよび Folded Gaussian
%  分布を構築する
%
% 入力:
%   ck_wnoise   - ノイズ付き信号系列 (1×N)
%   PAM_signal  - 送信シンボル系列 (1×N)\%   first_C     - first_C 値
%   first_F     - first_F 値
%   M           - モジュロ幅
%
% 出力 (構造体配列 fitResults):
%   fitResults(i).symbol  : 対象シンボル値
%   fitResults(i).mu      : 平均
%   fitResults(i).sigma   : 標準偏差
%   fitResults(i).xFold   : 折り返し分布の x 座標
%   fitResults(i).yFold   : 折り返し分布の PDF 値

% ck_wnoise_wTHP = ck_wnoise;
% M = bmax_value;
%%
    % トランケート
    idxStart = abs(first_C - first_F) + 1;
    ck_trunc = ck_wnoise(idxStart:end);
    PAM_trunc = PAM_signal(1:length(ck_trunc));
    % Mk

    symbols = unique(PAM_signal);
    P = 2*M;
    xFold = linspace(-M, M, 1000);

    fitResults = struct('symbol', [], 'mu', [], 'sigma', [], 'xFold', [], 'yFold', []);
    cnt = 0;
    for sym = symbols(:)'
        data = ck_trunc(PAM_trunc == sym) + 2*Mk(PAM_trunc == sym)*M; %M=bmax
        if isempty(data)
            continue;
        end
        cnt = cnt + 1;
        mu    = mean(data);
        sigma = std(data);
        % Folded Gaussian
        % ±4σ 範囲のシフト数
        nmin = ceil((mu - 4*sigma - M) / P);
        nmax = floor((mu + 4*sigma + M) / P);
        yFold = zeros(size(xFold));
        for n = nmin:nmax
            yFold = yFold + (1/(sigma*sqrt(2*pi))) * ...
                exp(-((xFold + n*P - mu).^2)/(2*sigma^2));
        end

        % fi=figure;
        % plot(xFold,yFold)
        % hold on
        % histogram(data,'Normalization','pdf')
        % close(fi) %絶対必要

        % 出力格納
        fitResults(cnt).symbol = sym;
        fitResults(cnt).mu     = mu;
        fitResults(cnt).sigma  = sigma;
        fitResults(cnt).xFold  = xFold;
        fitResults(cnt).yFold  = yFold;
    end
end
