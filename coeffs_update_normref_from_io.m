function [coeffs_eval, info] = coeffs_update_normref_from_io(coeffs_base, IO_run, tx_sym_run, win, tag)
% coeffs_update_normref_from_io
%   推奨案A: その run の送信列(tx_sym_run) と IO_run から LS を回して NormRef を推定し、
%   coeffs_base の NormRef だけを差し替えた coeffs_eval を返す。
%
% 入力:
%   coeffs_base : ふだん使っている係数（THP設計で使ったもの等）
%   IO_run      : その run の IO（y_wf, NoSpS, c_index を含む）
%   tx_sym_run  : その run で実際に送ったシンボル列（woTHPなら x_train、wTHPなら TX.bk）
%   win         : first_C, C_max, first_F, F_max
%   tag         : 表示用ラベル（省略可）
%
% 出力:
%   coeffs_eval : coeffs_base に対し NormRef のみ更新したもの
%   info        : NormRef_base, NormRef_run, gain_ratio など

    if nargin < 5
        tag = "";
    end

    % run固有の係数をLS推定（ここから NormRef を取り出す）
    coeffs_run = estimate_channel_ls(IO_run, tx_sym_run, win);

    coeffs_eval = coeffs_base;
    coeffs_eval.NormRef = coeffs_run.NormRef;

    info = struct();
    info.NormRef_base = coeffs_base.NormRef;
    info.NormRef_run  = coeffs_run.NormRef;

    if isfield(coeffs_base,'NormRef') && ~isempty(coeffs_base.NormRef) && abs(coeffs_base.NormRef) > 0
        info.gain_ratio = coeffs_run.NormRef / coeffs_base.NormRef;
    else
        info.gain_ratio = NaN;
    end

    % ログ（必要なら）
    if strlength(string(tag)) > 0
        gr = info.gain_ratio;
        fprintf('[%s] NormRef update: base=(%.4g%+.4gj), run=(%.4g%+.4gj), ratio=(%.4g%+.4gj), |ratio|=%.4g\n', ...
            string(tag), real(info.NormRef_base), imag(info.NormRef_base), ...
            real(info.NormRef_run),  imag(info.NormRef_run), ...
            real(gr), imag(gr), abs(gr));
    end
end
