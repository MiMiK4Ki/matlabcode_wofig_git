function saveCurrentFigure(f, foldername, figname)
    if ~exist(foldername, 'dir')
        mkdir(foldername);
    end
    
    base_path = fullfile(foldername, figname);
    % % pwd
    % base_path
    maxAttempts = 15;
    
    % savefig
    attempt = 0;
    while attempt < maxAttempts
        try
            savefig(f, base_path + ".fig")
            break;
        catch ME
            attempt = attempt + 1;
            fprintf('[savefig] %d回目の試行でエラー発生: %s\n', attempt, ME.message);
            if attempt == maxAttempts
                fprintf('[savefig] 最大試行回数に達したため中断します。\n');
                rethrow(ME);
            end
            pause(2);
        end
    end

    % exportgraphics (png)
    attempt = 0;
    while attempt < maxAttempts
        try
            exportgraphics(f, base_path + ".png");
            break;
        catch ME
            attempt = attempt + 1;
            fprintf('[exportgraphics_png] %d回目の試行でエラー発生: %s\n', attempt, ME.message);
            if attempt == maxAttempts
                fprintf('[exportgraphics_png] 最大試行回数に達したため中断します。\n');
                rethrow(ME);
            end
            pause(2);
        end
    end

    % exportgraphics (eps)
    attempt = 0;
    while attempt < maxAttempts
        try
            exportgraphics(f, base_path + ".eps");
            break;
        catch ME
            attempt = attempt + 1;
            fprintf('[exportgraphics_eps] %d回目の試行でエラー発生: %s\n', attempt, ME.message);
            if attempt == maxAttempts
                fprintf('[exportgraphics_eps] 最大試行回数に達したため中断します。\n');
                rethrow(ME);
            end
            pause(2);
        end
    end
    
    close(f)
end
