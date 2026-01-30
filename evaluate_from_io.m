function [Res, D] = evaluate_from_io(IO, coeffs, param, ref_sym, mode)
% evaluate_from_io
%   連続波形 IO.y_wf をサンプリングして y_sym を作り、
%   evaluate_from_sym に渡して評価する。
%
% 必須:
%   IO.y_wf, IO.NoSpS, IO.c_index
%   coeffs.first_C, coeffs.NormRef

y  = IO.y_wf(:).';
L  = IO.NoSpS;
c_index = IO.c_index;

x_ref = ref_sym(:).';

% 取り出せるシンボル数（c_index基準）
NsymAvail = min(numel(x_ref), floor((numel(y)-c_index)/L) + 1);

% 旧main準拠: first_C を含めた位置でサンプルを切る
idx = c_index + (coeffs.first_C + (0:NsymAvail-1))*L;

% 離散系列（NormRefで正規化）
y_sym_raw = y(idx);
y_sym     = y_sym_raw ./ coeffs.NormRef;

% 評価（EVM含む）
[Res, Dsym] = evaluate_from_sym(y_sym, coeffs, param, x_ref(1:NsymAvail), mode);

% debug
if nargout >= 2
    D = Dsym;
    D.idx       = idx;
    D.y_sym     = y_sym;
    D.y_sym_raw = y_sym_raw;
    D.c_index   = c_index;
    D.NoSpS     = L;
else
    D = [];
end

end
