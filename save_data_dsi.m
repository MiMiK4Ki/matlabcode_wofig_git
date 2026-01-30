function save_data_dsi(foldername, first_F_value, suffix, data)
    % 保存するファイルパスの作成
    filepath = fullfile(foldername, strcat("vardsi", suffix, "fF", num2str(first_F_value), ".mat"));
    
    save(filepath, '-struct', 'data');
    

end
