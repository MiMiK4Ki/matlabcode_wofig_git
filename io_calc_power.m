function [P_ref, info] = io_calc_power(IO, K, param, method)
% io_calc_power  参照電力 P_ref を計算
%
% method:
%   "continuous" : 連続波形 y_wf を区間積分（旧 calculate_Power と同じ系）
%   "zoh"        : 記号中心サンプルを矩形ホールドしたとみなす近似
%
% 使い方例:
%   K = numel(x_sym);
%   [P_ref,info] = io_calc_power(IO, K, param, "continuous");

if nargin < 4 || isempty(method)
    method = "continuous";
end

y = IO.y_wf(:).';
L = IO.NoSpS;
c = IO.c_index;

% Tc と dt（実験なら IO.Tc / IO.fs を優先）
if isfield(IO,'Tc') && ~isempty(IO.Tc)
    Tc = IO.Tc;
else
    Tc = param.Tc_but / param.cutoff_coeff;
end

if isfield(IO,'fs') && ~isempty(IO.fs)
    dt = 1 / IO.fs;
else
    dt = Tc / L;
end

if method == "continuous"
    s1 = c;
    s2 = c + (K-1)*L;
    seg = y(s1:s2);

    % 連続平均電力（複素でもOK）
    P_ref = sum(abs(seg).^2) * dt / ((K-1)*Tc);

    info = struct('method',"continuous",'dt',dt,'Tc',Tc,'s1',s1,'s2',s2,'K',K);

elseif method == "zoh"
    % 記号中心サンプル列（この列を Tc 区間ホールドしたとみなす）
    idx = c + (0:K-1)*L;
    y_sym = y(idx);

    % 矩形近似なら平均電力は mean(|y_sym|^2) と一致
    P_ref = mean(abs(y_sym).^2);

    info = struct('method',"zoh",'dt',dt,'Tc',Tc,'idx',idx,'K',K);

else
    % 余計なエラー処理は増やさない方針なので、ここは最小限
    error('io_calc_power: method must be "continuous" or "zoh".');
end

end
