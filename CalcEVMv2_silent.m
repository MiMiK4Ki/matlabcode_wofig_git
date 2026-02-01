function [rx_fit, err, EVM, EVMref, info] = CalcEVMv2_silent(rx, tx, M)
% CalcEVMv2_silent
%   PAM/QAM 両対応の EVM(v2) を返す（表示なし）
%
% 目的（意図）:
%   1) 受信点 rx を 1係数の複素ゲインで tx に最小二乗フィット（残留の回転/ゲインを除去）
%   2) 誤差 err = rx_fit - tx の RMS を計算
%   3) EVM:    tx の平均電力で正規化
%      EVMref: 理想星座（PAM/QAM）の平均電力で正規化
%
% 出力:
%   rx_fit : best-fit 後の受信点
%   err    : error vector
%   EVM    : [%]
%   EVMref : [%]（理想星座電力正規化）
%   info   : 付帯情報（fit係数など）

    if nargin < 3
        M = [];
    end

    rx = rx(:);
    tx = tx(:);

    N = min(numel(rx), numel(tx));
    rx = rx(1:N);
    tx = tx(1:N);

    % NaN/Inf 除去
    ok = isfinite(real(rx)) & isfinite(imag(rx)) & isfinite(real(tx)) & isfinite(imag(tx));
    rx = rx(ok);
    tx = tx(ok);

    if isempty(rx)
        rx_fit = rx; err = rx;
        EVM = NaN; EVMref = NaN;
        info = struct('alpha',NaN,'Ptx',NaN,'Pconst',NaN,'N',0);
        return;
    end

    % --- best-fit 1係数（残留の定数複素ゲインを除去） ---
    % 目的: alpha*rx ≈ tx を最小二乗で満たす alpha
    denom = (rx' * rx);
    if abs(denom) < eps
        alpha = 1;
    else
        alpha = (rx' * tx) / denom;   % alpha = (rx^H tx)/(rx^H rx)
    end

    rx_fit = alpha * rx;
    err = rx_fit - tx;

    % --- RMS ---
    Pe = mean(abs(err).^2);
    Ptx = mean(abs(tx).^2);

    if Ptx <= 0
        EVM = NaN;
    else
        EVM = 100 * sqrt(Pe / Ptx);
    end

    % --- EVMref: 理想星座の平均電力で正規化 ---
    Pconst = local_constellation_avg_power(tx, M);
    if isempty(Pconst) || ~isfinite(Pconst) || Pconst <= 0
        EVMref = NaN;
    else
        EVMref = 100 * sqrt(Pe / Pconst);
    end

    info = struct();
    info.alpha = alpha;
    info.Pe    = Pe;
    info.Ptx   = Ptx;
    info.Pconst = Pconst;
    info.N = numel(tx);
end

% ===== local helper =====

function Pconst = local_constellation_avg_power(tx, M)
% tx が実数なら PAM、複素なら QAM とみなす（PAMでも rx が複素になり得るが tx は実数のはず）
% 可能なら tx に出現するユニーク点から電力を出す（スケールが変わっても追従できる）
    if nargin < 2 || isempty(M)
        M = numel(unique(tx));
    end

    if ~isscalar(M) || M < 2 || ~isfinite(M)
        Pconst = NaN;
        return;
    end

    ux = unique(tx);
    % 十分な長さの系列なら通常 ux は全点を含むので、それを理想星座とみなす
    if numel(ux) == M
        Pconst = mean(abs(ux).^2);
        return;
    end

    if isreal(tx)
        % PAM（標準: levels=-(M-1):2:(M-1) を想定）
        levels = (-(M-1):2:(M-1));
        Pconst = mean(levels.^2);
        return;
    else
        % QAM（正方形）
        L = sqrt(double(M));
        if abs(L - round(L)) > 1e-12
            Pconst = NaN;
            return;
        end
        L = round(L);
        levels = (-(L-1):2:(L-1));
        [I,Q] = meshgrid(levels, levels);
        const = I(:) + 1j*Q(:);
        Pconst = mean(abs(const).^2);
        return;
    end
end
