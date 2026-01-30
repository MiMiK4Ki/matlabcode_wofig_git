function [tap, getTap] = make_tap(enabled)
% 使い方:
%   [tap, getTap] = make_tap(true);
%   tap('foo', 123); S = getTap();

    if ~enabled
        tap    = @(varargin) [];     % 何もしない
        getTap = @() struct();       % 空を返す
        return
    end

    S = struct();                    % ここに貯める

    % --- "記録する" 側（S に追記） ---
    function capture(name, val)
        fld = matlab.lang.makeValidName(char(name));
        S.(fld) = shrink(val);
    end

    % --- "取り出す" 側（最新の S を返す） ---
    function Sout = pull()
        Sout = S;
    end

    tap    = @capture;
    getTap = @pull;

    function v = shrink(v)
        N = 2000;                    % 研究用の軽い縮約（任意）
        if isnumeric(v) && numel(v) > N, v = v(1:N); end
    end
end
