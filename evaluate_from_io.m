function [Res, D] = evaluate_from_io(IO, coeffs, param, ref_sym, mode)
% evaluate_from_io
%   連続波形 IO.y_wf をサンプリングして y_sym を作り、
%   evaluate_from_sym に渡して評価する。
%
% 必須:
%   IO.y_wf, IO.NoSpS, IO.c_index
%   coeffs.first_C, coeffs.NormRef

% 連続波形 → 離散シンボル
Ddisc = io_discretize(IO, coeffs, ref_sym);
y_sym_raw = Ddisc.y_sym(:).';
y_sym     = y_sym_raw ./ coeffs.NormRef;

% 評価（EVM含む）
[Res, Dsym] = evaluate_from_sym(y_sym, coeffs, param, ref_sym(1:Ddisc.K), mode);

% debug
if nargout >= 2
    D = Dsym;
    D.idx       = Ddisc.idx;
    D.y_sym     = y_sym;
    D.y_sym_raw = y_sym_raw;
    D.c_index   = Ddisc.c_index;
    D.NoSpS     = IO.NoSpS;
else
    D = [];
end

end
