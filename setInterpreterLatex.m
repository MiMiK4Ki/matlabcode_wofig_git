function setInterpreterLatex(fontSize)
    % デフォルトのフォントサイズを設定
    if nargin < 1
        fontSize = 14;
    end

    grid on
    
    % 現在の軸を取得
    ax = gca;
    
    % インタープリターをlatexに設定
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.Title.Interpreter = 'latex'; % タイトルがある場合
    if ~isempty(ax.Legend) % 凡例がある場合
        ax.Legend.Interpreter = 'latex';
    end
    
    % フォントサイズを設定
    ax.FontSize = fontSize;
end
