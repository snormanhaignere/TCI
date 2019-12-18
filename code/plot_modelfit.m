function plot_modelfit(diff_context, same_context, diff_context_bestpred, ...
    loss_best_model, best_intper_sec, best_delay_sec_median, best_shape, fname_global, ...
    intper_sec, delay_sec_start, unique_segs, lag_t,...
    plot_win, plot_smoothwin, plot_delaystat, plot_delay_range, ploterrquant, linewidth, figh)

corr_range = quantile(diff_context(:), [0.01, 0.99]);

win_string = sprintf('-win-%d-%dms', round(plot_win(1)*1000), round(plot_win(2)*1000));

if plot_smoothwin>0    
    same_context = mysmooth(same_context, 1/diff(lag_t(1:2)), plot_smoothwin*1000);
    diff_context = mysmooth(diff_context, 1/diff(lag_t(1:2)), plot_smoothwin*1000);
    diff_context_bestpred = mysmooth(diff_context_bestpred, 1/diff(lag_t(1:2)), plot_smoothwin*1000);
    win_string = [win_string sprintf('-smooth-%dms', round(plot_smoothwin*1000))];
end

% plot prediction for best delay, lineplot
if ~isempty(diff_context_bestpred)
    
    clf(figh);
    set(figh, 'Position', [100, 100, 900, 900]);
    
    clear X;
    invariance_line = NaN;
    for k = 1:size(diff_context,2)
        subplot(4, 2, k);
        hold on;
        plot(lag_t([1 end]) * 1000, [0,0], 'k--', 'LineWidth', linewidth);
        h1 = plot(lag_t * 1000, same_context(:,k), 'LineWidth', linewidth);
        h2 = plot(lag_t * 1000, diff_context(:,k), 'LineWidth', linewidth);
        h3 = plot(lag_t * 1000, diff_context_bestpred(:,k), 'LineWidth', linewidth);
        plot(unique_segs(k)*[1 1], corr_range, 'k--', 'LineWidth', linewidth);
        if ~isnan(invariance_line); plot(plot_win * 1000, invariance_line*[1 1], 'k--', 'LineWidth', 2); end
        xlim(plot_win * 1000);
        ylim(corr_range);
        xlabel('Time Lag (ms)');
        ylabel('Pearson Correlation');
        title(sprintf('Seg: %.0f ms', unique_segs(k)))
        if k == 1
            legend([h1, h2, h3], 'Same', 'Cross', 'Model');
        end
        box off;
    end
    fname = [fname_global '-prediction-lineplot' win_string];
    export_fig([fname '.pdf'], '-pdf', '-transparent');
    export_fig([fname '.png'], '-png', '-transparent', '-r150');
    savefig(figh, [fname '.fig']);
    
    % plot prediction for best delay, image
    clf(figh);
    set(figh, 'Position', [100, 100, 900, 600]);
    subplot(2,1,1);
    imagesc(diff_context', corr_range(2) * [-1, 1]);
    subplot(2,1,2);
    imagesc(diff_context_bestpred', corr_range(2) * [-1, 1]);
    for i = 1:2
        subplot(2,1,i);
        colorbar;
        colormap(flipud(cbrewer('div', 'RdBu', 128)));
        set(gca, 'YTick', 1:length(unique_segs), 'YTickLabel', unique_segs);
        xticks = get(gca, 'XTick');
        set(gca, 'XTick', xticks, 'XTickLabel', lag_t(xticks)*1000);
        ylabel('Seg Dur (ms)'); xlabel('Lag (ms)');
        if i == 2
            title(sprintf('rf=%.f ms, delay=%.f ms', best_intper_sec*1000, best_delay_sec_median*1000));
        end
    end
    fname = [fname_global '-prediction-image' win_string];
    export_fig([fname '.png'], '-png', '-transparent', '-r150');
    savefig(figh, [fname '.fig']);
    
end

% delays to plot
X_start_delays = loss_best_model;
if ~strcmp(plot_delaystat, 'start')
    delay_sec_altstat = nan(size(X_start_delays));
    for l = 1:length(intper_sec)
        delay_sec_altstat(l,:) = modelwin_convert_delay(intper_sec(l), delay_sec_start, best_shape, plot_delaystat);
    end
    delays_to_plot = plot_delay_range(1):diff(delay_sec_start(1:2)):plot_delay_range(2);
    X_altdelays = nan(length(intper_sec), length(delays_to_plot));
    for l = 1:length(intper_sec)
        xi = delays_to_plot > delay_sec_altstat(l,1) & delays_to_plot < delay_sec_altstat(l,end);
        X_altdelays(l,xi) = interp1(delay_sec_altstat(l,:), X_start_delays(l,:), delays_to_plot(xi));
    end
    X_altdelays(isnan(X_altdelays)) = inf;
    X = X_altdelays;
else
    X = X_start_delays;
end

% plot the error vs. parameters
cmap = flipud(cbrewer('seq', 'Reds', 128));
[minX,zi] = min(X(:));
[~, xi] = ind2sub(size(X), zi);
bounds = [minX, quantile(X(:,xi), ploterrquant)];
clear xi zi;
if ~all(isnan(X(:)))
    clf(figh);
    set(figh, 'Position', [100, 100, 600, 600]);
    imagesc(X, bounds);
    colormap(cmap);
    colorbar;
    yticks = unique(round(linspace(1, length(intper_sec), 5)));
    set(gca, 'YTick', yticks, 'YTickLabel', 1000*intper_sec(yticks));
    xticks = unique(round(linspace(1, length(delays_to_plot), 5)));
    set(gca, 'XTick', xticks, 'XTickLabel', 1000*delays_to_plot(xticks));
    xlabel('Delay (ms)'); ylabel('Receptive Field Duration (ms)');
    title(sprintf('rf=%.f ms, delay=%.f ms', best_intper_sec*1000, best_delay_sec_median*1000));
    set(gca, 'FontSize', 12);
    fname = [fname_global '-model-error'];
    export_fig([fname '.png'], '-png', '-transparent', '-r150');
    savefig(figh, [fname '.fig']);
end