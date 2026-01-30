function plot_io_overlay(IO, ChanFrame, K_plot)
% 連続波形とシンボル中心サンプル、スケール済み x_sym を1枚に重ねて可視化
% K_plot: 図示する先頭シンボル数（省略時 200）

    if nargin < 3 || isempty(K_plot), K_plot = min(200, numel(IO.x_sym)); end

    dt    = ChanFrame.dt;
    NoSpS = IO.NoSpS;
    y     = IO.y_wf(:).';
    x     = IO.x_sym(:).';

    % 主タップ位置（最大振幅インデックス）
    [c, c_index] = max(abs(ChanFrame.h_rx));

    % シンボル中心サンプルの抽出インデックス
    idx = c_index + (0:K_plot-1)*NoSpS;
    idx = idx(idx <= numel(y));
    t   = (0:numel(y)-1) * dt;
    tc  = (idx-1) * dt;     % center sample の時刻

    % 表示用ゲイン（最初の非ゼロシンボルから推定）
    gain =c;

    figure; hold on;
    plot(t, y, 'DisplayName','y\_wf (continuous)');
    plot(tc, y(idx), 'o', 'DisplayName','y @ symbol centers');
    stem(tc(1:numel(idx)), x(1:numel(idx))*gain, 'filled', 'DisplayName','x\_sym (scaled)');
    grid on; xlabel('Time [s]'); ylabel('Amplitude');
    title(sprintf('Overlay: continuous y and symbol centers (K\\_plot=%d)', numel(idx)));
    legend('show','Location','best');
end
