function parSave(filename, iterationRes)
    % 外部関数でsaveを呼び出すことで、parfor内での直接呼び出しの制限を回避します。
    save(filename, 'iterationRes', '-v7.3');
end