function ChanFrame = chan_coeffs_to_chanframe_interp(param, coeffs, win, method)
% chan_coeffs_to_chanframe_interp
%   推定した離散 ChannelCoeff から「近似の連続IR」を作る（io_make_theory 用）
%   ※厳密な連続IR復元ではなく、シンボル間を補間して埋める簡易版
%
% 入力:
%   param.NoSpS, param.Tc_but, param.cutoff_coeff
%   coeffs.ChannelCoeff, coeffs.NormRef
%   win.first_C, win.C_max
%   method: "linear" or "previous" etc (interp1 の方法)
%
% 出力:
%   ChanFrame.h_rx, ChanFrame.dt

if nargin < 4
    method = "linear";
end

NoSpS = param.NoSpS;
Tc = param.Tc_but / param.cutoff_coeff;
dt = Tc / NoSpS;

m = win.first_C:win.C_max;

% ChannelCoeff は NormRef で正規化されている前提 → 元スケールに戻す
h_sym = (coeffs.ChannelCoeff(:).') * coeffs.NormRef;

% シンボル位置のサンプルindex（m=first_C → 1）
idxGrid = 1 + (m - win.first_C)*NoSpS;

% 補間で連続サンプル列を生成
n = 1:idxGrid(end);
h_rx = interp1(idxGrid, h_sym, n, method);

% 右側を少しだけゼロパディング（畳み込みの尻尾用）
h_rx = [h_rx, zeros(1, 5*NoSpS)];

ChanFrame = struct();
ChanFrame.h_rx = h_rx;
ChanFrame.dt   = dt;
end
