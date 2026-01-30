function symbols = generate_symbols(param, TXD_N)
% generate_symbols
%   param.MODTYPE に応じて PAM/QAM シンボル列を生成する。
%   既定: MODTYPE が未指定なら "PAM" 扱い。

    if nargin < 2 || isempty(TXD_N)
        TXD_N = param.TXD_N;
    end

    modtype = "PAM";
    if isfield(param, "MODTYPE") && ~isempty(param.MODTYPE)
        modtype = upper(string(param.MODTYPE));
    end

    switch modtype
        case "QAM"
            symbols = generate_QAM(param.MODNUM, TXD_N);
        otherwise
            symbols = generate_PAM(param.MODNUM, TXD_N);
    end
end
