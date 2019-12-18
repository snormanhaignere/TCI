function corr_range = plot_cross_context_corr(...
    same_context, diff_context, unique_segs, lag_t, ...
    plot_range, plot_win, simfunc, figh)

% plot
clf(figh);
set(figh, 'Position', [100, 100, 900, 900]);
if isnan(plot_range)
    ti = lag_t >= plot_win(1) & lag_t <= plot_win(2);
    X = cat(3, same_context(ti,:), diff_context(ti,:));
    if any(ismember({'mae'}, simfunc))
        corr_range = quantile(X(X~=0), [0.01, 0.99]);
    else
        corr_range = [-1 1] * max(X(:))*1.05;
    end
    clear ti X;
else
    corr_range = plot_range;
end
n_seg_durs = size(same_context,2);
for k = 1:n_seg_durs
    subplot(4, 2, k);
    X = [same_context(:,k), diff_context(:,k)];
    plot(lag_t * 1000, X, 'LineWidth', 2);
    hold on;
    plot([0 0], corr_range, 'r--', 'LineWidth', 2);
    plot(unique_segs(k)*[1 1], corr_range, 'r--', 'LineWidth', 2);
    plot(plot_win * 1000, [0 0], 'k--', 'LineWidth', 2);
    xlim(plot_win * 1000);
    ylim(corr_range);
    xlabel('Time Lag (ms)');
    ylabel('Similarity');
    title(sprintf('Seg: %.0f ms', unique_segs(k)))
    if k == 1
        legend('Same', 'Cross');
    end
end