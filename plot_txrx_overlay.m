function plot_txrx_overlay(IO, coeffs, param, x_train, mode, tag, K_plot, tx_dk, lin_conv_output, K_hist)
% plot_txrx_overlay
% 送信シンボル・受信シンボル結果・連続受信波形（正規化）を重ね描きして、
% サンプリングタイミングとズレ(first_C/first_F)の整合を目視確認する。
%
% ★追加機能:
%   - 送信シンボル値 x_train ごとの RX_soft ヒストグラム（pdf）を表示
%   - K_plot とは別に K_hist を用意（デフォルトで全体 ~20000点）
%
% 入力:
%   IO      : struct, 必須 .y_wf .NoSpS .c_index
%   coeffs  : struct, 必須 .NormRef .first_C .first_F .ChannelCoeff
%   param   : struct, 必須 .MODNUM .bmax_initial
%   x_train : 送信シンボル列（参照）
%   mode    : 'woTHP' or 'wTHP'
%   tag     : 図タイトル用文字列（例 'A5 Theory woTHP'）
%   K_plot  : 表示する「比較対象TXシンボル数」（デフォルト2000）
%   tx_dk   : (wTHP用) 追加で重ねたい送信列（例 TX.dk）。woTHPなら省略/[]でOK
%   lin_conv_output : (任意) conv(x,h) などの検証用波形（idx_common に対応する長さ推奨）
%   K_hist  : (任意) ヒストグラム用のTXシンボル数（デフォルト ~20000）
%
% 備考:
%   - RX_soft ヒストグラムは、(TXと整列した) rx_soft_common を使用する
%   - pdf 表示のため histogram(...,'Normalization','pdf') を用いる

    if nargin < 6 || isempty(tag),    tag = ""; end
    if nargin < 7 || isempty(K_plot), K_plot = min(2000, numel(x_train)); end
    if nargin < 8, tx_dk = []; end
    if nargin < 9, lin_conv_output = []; end
    if nargin < 10, K_hist = []; end

    % ---- 互換: 9番目がスカラーで10番目が無い場合は K_hist とみなす ----
    if isempty(K_hist)
        if ~isempty(lin_conv_output) && isscalar(lin_conv_output)
            K_hist = lin_conv_output;
            lin_conv_output = [];
        else
            K_hist = min(20000, numel(x_train));
        end
    else
        if isempty(K_hist) || ~isfinite(K_hist) || K_hist <= 0
            K_hist = min(20000, numel(x_train));
        end
    end

    mode = string(mode);

    y       = IO.y_wf(:).';
    NoSpS   = IO.NoSpS;
    c_index = IO.c_index;

    first_C = coeffs.first_C;
    first_F = coeffs.first_F;
    NormRef = coeffs.NormRef;

    % ===== 連続波形の正規化 =====
    idx_m0 = 1 - first_C;                 % m=0 のインデックス
    A0     = coeffs.ChannelCoeff(idx_m0); % m=0 の係数

    % ===== シフト量（旧規約に準拠） & スケール =====
    if mode == "woTHP"
        shift_sym = abs(first_C);
        scale = (NormRef * A0);
    else
        shift_sym = abs(first_C - first_F);
        scale = NormRef;
    end

    yN = y ./ scale;

    % ===== 受信側サンプル列の長さ上限（evaluate_from_io と同型）=====
    Krx_max = min(numel(x_train), floor((numel(yN) - c_index) / NoSpS) + 1);

    % ============================================================
    % (1) オーバレイ表示用（K_plot）
    % ============================================================
    K_show = min(K_plot, numel(x_train));
    K_use  = min(Krx_max, K_show + shift_sym);

    idx_rx = c_index + (first_C + (0:K_use-1)) * NoSpS;
    y_smp  = yN(idx_rx);

    Mmod       = param.MODNUM;
    bmax_value = param.bmax_initial * (Mmod/2);

    switch mode
        case "woTHP"
            rx_soft = y_smp;
        case "wTHP"
            rx_soft = wrapToM(y_smp, bmax_value);
        otherwise
            error("mode must be 'woTHP' or 'wTHP'.");
    end

    modtype = "PAM";
    if isfield(param, "MODTYPE") && ~isempty(param.MODTYPE)
        modtype = upper(string(param.MODTYPE));
    elseif ~isreal(x_train)
        modtype = "QAM";
    end

    % hard decision（PAM/QAM）
    switch modtype
        case "QAM"
            [RX_sym_idx, const] = qam_to_symbols(rx_soft, Mmod);
            rx_hard = const(RX_sym_idx);
        otherwise
            [RX_sym_idx, ~] = pam_to_symbols(real(rx_soft), Mmod);
            rx_hard = levels_from_sym(RX_sym_idx, Mmod);
    end

    % TX と RX を同じ横軸に合わせる
    Kc = K_use - shift_sym;

    if ~isempty(tx_dk)
        tx_dk = tx_dk(:).';
        Kc = min(Kc, numel(tx_dk));
    end

    idx_common      = idx_rx(shift_sym+1 : shift_sym+Kc);
    tx_train_common = x_train(1:Kc);

    y_smp_common    = y_smp(shift_sym+1 : shift_sym+Kc);
    rx_soft_common  = rx_soft(shift_sym+1 : shift_sym+Kc);
    rx_hard_common  = rx_hard(shift_sym+1 : shift_sym+Kc);

    if ~isempty(tx_dk)
        tx_dk_common = tx_dk(1:Kc);
    else
        tx_dk_common = [];
    end

    % 連続波形の表示範囲（見やすさ優先）
    n1 = max(1, idx_rx(1) - 2*NoSpS);
    n2 = min(numel(yN), idx_rx(end) + 2*NoSpS);
    n  = n1:n2;

    % ===== オーバレイプロット（実部/虚部で分離）=====
    local_plot_component("Real", @real, n, yN, idx_common, y_smp_common, rx_soft_common, rx_hard_common, ...
        tx_train_common, tx_dk_common, lin_conv_output, tag, mode, Kc, shift_sym);
    if ~isreal(yN) || ~isreal(tx_train_common) || ~isreal(rx_soft_common)
        local_plot_component("Imag", @imag, n, yN, idx_common, y_smp_common, rx_soft_common, rx_hard_common, ...
            tx_train_common, tx_dk_common, lin_conv_output, tag, mode, Kc, shift_sym);
    end

    % ============================================================
    % (2) ヒストグラム用（K_hist）：TX値ごとの RX_soft 分布（pdf）
    % ============================================================
    K_hist = min(round(K_hist), numel(x_train));
    K_use_h = min(Krx_max, K_hist + shift_sym);

    if K_use_h <= shift_sym
        warning('plot_txrx_overlay: K_hist too small after shift. Skip histogram.');
        return;
    end

    idx_rx_h = c_index + (first_C + (0:K_use_h-1)) * NoSpS;
    y_smp_h  = yN(idx_rx_h);

    switch mode
        case "woTHP"
            rx_soft_h = y_smp_h;
        case "wTHP"
            rx_soft_h = wrapToM(y_smp_h, bmax_value);
    end

    Kc_h = K_use_h - shift_sym;
    local_plot_hist("Real", real(x_train(1:Kc_h)).', real(rx_soft_h(shift_sym+1 : shift_sym+Kc_h)).', tag, mode);
    if ~isreal(rx_soft_h) || ~isreal(x_train)
        local_plot_hist("Imag", imag(x_train(1:Kc_h)).', imag(rx_soft_h(shift_sym+1 : shift_sym+Kc_h)).', tag, mode);
    end

end


function lv = levels_from_sym(sym_idx, Mmod)
    lv_ideal = (-(Mmod-1):2:(Mmod-1));
    lv = lv_ideal(sym_idx);
end

function local_plot_component(tag_comp, comp, n, yN, idx_common, y_smp_common, rx_soft_common, rx_hard_common, ...
    tx_train_common, tx_dk_common, lin_conv_output, tag, mode, Kc, shift_sym)
    figure; hold on;

    plot(n, comp(yN(n)), 'DisplayName',sprintf('y\\_wf / scale (%s)', tag_comp));
    plot(idx_common, comp(y_smp_common), '.', 'DisplayName',sprintf('samples (y\\_smp, %s)', tag_comp));
    plot(idx_common, comp(rx_soft_common), 'x', 'DisplayName',sprintf('RX soft (%s)', tag_comp));
    stem(idx_common, comp(rx_hard_common), 'filled', 'DisplayName',sprintf('RX hard (%s)', tag_comp));
    stem(idx_common, comp(tx_train_common), 'DisplayName',sprintf('TX x\\_train (%s)', tag_comp));

    if ~isempty(tx_dk_common)
        stem(idx_common, comp(tx_dk_common), 'DisplayName',sprintf('TX dk (%s)', tag_comp));
    end

    if ~isempty(lin_conv_output)
        stem(idx_common, comp(lin_conv_output(1:numel(idx_common))), 'DisplayName',sprintf('conv(x,h) (%s)', tag_comp));
    end

    grid on;
    xlabel('sample index');
    ylabel('amplitude (normalized)');
    title(sprintf('%s  [%s]  %s (K=%d, shift=%d)', tag, mode, tag_comp, Kc, shift_sym));
    legend('show','Location','best');
end

function local_plot_hist(tag_comp, tx_h, rx_h, tag, mode)
    mu = mean(rx_h);
    sg = std(rx_h);
    if ~isfinite(sg) || (sg < 1e-12)
        sg = max(abs(rx_h - mu)) + 1e-12;
    end
    rmin = mu - 6*sg;
    rmax = mu + 6*sg;
    if rmin == rmax
        rmin = mu - 1;
        rmax = mu + 1;
    end
    Nb = 80;
    edges = linspace(rmin, rmax, Nb+1);

    lv = unique(tx_h(:));
    lv = sort(lv);

    figure; hold on; grid on;
    for ii = 1:numel(lv)
        mask = (tx_h == lv(ii));
        data = rx_h(mask);
        if isempty(data), continue; end
        histogram(data, edges, ...
            'Normalization','pdf', ...
            'DisplayStyle','stairs', ...
            'LineWidth',1.5, ...
            'DisplayName',sprintf('Tx=%g (N=%d)', lv(ii), sum(mask)));
    end

    xlabel(sprintf('RX soft amplitude (%s)', tag_comp));
    ylabel('pdf');
    title(sprintf('%s  [%s]  RX soft pdf per Tx (%s)  (N=%d)', ...
        tag, mode, tag_comp, numel(rx_h)));
    legend('show','Location','best');

    e = rx_h - tx_h;
    figure;
    histogram(e, 'Normalization','pdf');
    grid on;
    xlabel(sprintf('RX\\_soft - TX (%s)', tag_comp));
    ylabel('pdf');
    title(sprintf('%s  [%s]  Error pdf (%s, var=%.6g)  (N=%d)', ...
        tag, mode, tag_comp, var(e), numel(e)));
end
