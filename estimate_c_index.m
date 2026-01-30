function [c_index, info] = estimate_c_index(IO, x_sym, varargin)
% estimate_c_index  x と連続 y だけで c_index を推定（FIR を同時に仮推定）
% 入力:
%   IO     : struct('y_wf','NoSpS')   % 収録波形とオーバーサンプリング比
%   x_sym  : 送信記号列（M-PAM レベル）
% オプション:
%   'Q'        : FIR次数（既定 15）
%   'Smax'     : 記号ずれ探索範囲（既定 200）
%   'skipSyms' : 先頭捨てる記号数（既定 20; 立上り除去）
%
% 出力:
%   c_index : 1始まりの中心サンプル位置
%   info    : 位相/ずれ/誤差など

    p = inputParser;
    p.addParameter('Q',15,@(v)isnumeric(v)&&isscalar(v)&&v>=1);
    p.addParameter('Smax',200,@(v)isnumeric(v)&&isscalar(v)&&v>=0);
    p.addParameter('skipSyms',20,@(v)isnumeric(v)&&isscalar(v)&&v>=0);
    p.parse(varargin{:});
    Q        = p.Results.Q;
    Smax     = p.Results.Smax;
    skipSyms = p.Results.skipSyms;

    L  = IO.NoSpS;
    y  = IO.y_wf(:).';
    x  = x_sym(:).';

    best.err = Inf; best.p = 1; best.s = 0; best.Q = Q; best.Kuse = 0;

    for p0 = 1:L
        ys = y(p0:L:end);            % 位相 p0 のサブサンプリング列
        Lp = numel(ys);
        Smax_eff = min(Smax, max(0, Lp - (Q+50))); % 行数確保のための上限

        for s = 0:Smax_eff
            Y = ys(1+s:end);                     % 記号ずれ s を適用
            Kmax = min(numel(x), numel(Y));
            Keff = Kmax - (Q-1) - skipSyms;
            if Keff < 50, continue; end          % 最低点数の確保

            n = (skipSyms + Q) : (skipSyms + Q + Keff - 1);
            Yseg = Y(n);
            X = zeros(Keff, Q);
            for m = 0:Q-1
                X(:,m+1) = x(n - m);
            end

            % DC除去（オフセット抑制）
            Yz = Yseg - mean(Yseg);
            Xz = X - mean(X,1);

            % LS（必要なら正則化に差し替え可）
            h = Xz \ Yz;
            err = mean( (Yz - Xz*h).^2 );

            if err < best.err
                best.err = err; best.p = p0; best.s = s;
                best.h = h; best.Kuse = Keff;
            end
        end
    end

    c_index = best.p + best.s * L;
    info = struct('method','phase+LS','phase',best.p,'sym_shift',best.s, ...
                  'MSE',best.err,'Q',best.Q,'Kuse',best.Kuse);
end
