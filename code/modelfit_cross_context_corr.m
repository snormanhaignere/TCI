function [M, MAT_file] = modelfit_cross_context_corr(L, varargin)

% Fit parametric window using cross-context correlation data. This function
% should be applied to the output structure L returned by
% cross_context_corr.m

%% Default parameters

clear I;

% distribution used to model the window
% 'gamma' or 'gauss'
I.distr = 'gamma';

% shape parameter for parametric window, not relevant for gaussian
% 1 -> exponential, higher -> more gaussian
I.shape = 1;

% range of integration periods to consider in seconds
I.intper_range = L.unique_segs([1 end])/1000; % in seconds

% number of integration periods to consider between this range
% spaced on a logarithmic scale
I.nintper = 20;

% how to calculate delay
% options are 'peak', 'median', or 'start' (if causal)
I.delaypoint = 'median';

% range of delays to consider in seconds
I.delay_range = [0, L.unique_segs(end)/1000];

% interval between delays in seconds
I.delay_interval = 1/L.sr;

% whether and by how much to ramp segments
% this should match the actual ramps used 
% in the stitching together of segments
I.rampwin = L.rampwin;
I.rampdur = L.rampdur; % in seconds

% whether to normalize cross-context correlation
% divisively by same-context correlation for each lag
I.divnorm = true;

% lags with a same-context correlation below this value are excluded
% the same-context correlation which appears in the denominator 
% (i.e. cross-context / same-context)
I.mindenom = 0.01;

% an alternative to divcorr
% here we normalize by the average of the same-context
% correlation over the window
I.normcorr = false;

% whether to weight the errors by the strength
% of the same-context correlation which appears
% in the denomitor (i.e. cross-context / same-context)
I.weightdenom = true;

% range of shape parameters for gamma distribution
% 1 -> exponential, higher -> more Gaussian
I.shape = logspace(log10(1), log10(10), 10);

% whether to use the power ratio between 
% segments to predict the correlation
I.winpowratio = true;

% returns the best causal window, excluding non-causal windows
I.bestcausal = true;

% whether to force all windows to be causal
% by zeroing time-points before stim onset
I.forcecausal = false;

% whether to only show errors in plots for causal windows
I.plotcausal = true;

% whether to transform the weighting applied to segments
I.tranweightnsegs = 'none'; % applied to total segs

% whether to transform the weighting applied to the denominator
I.tranweightdenom = 'none';

% whether to transform the predictions
I.tranpred = 'none';

% whether to overwrite existing results
I.overwrite = false;

% whether to enter debug mode
I.keyboard = false;

% whether to plot figures
I.plot_figure = true;

% window used for plotting
I.plot_win = L.lag_t([1 end]); % in seconds

% whether to run the analysis
% if false, just returns the parameters to be used
I.run = true;

% figure handle used for plotting
I.figh = matlab.ui.Figure.empty;

%% Parse user-specified parameters, create parameter string for this analysis

[I, ~, C_value] = parse_optInputs_keyvalue(varargin, I);

% if using gaussian distribution
% shape is fixed
if strcmp(I.distr, 'gauss')
    I.shape = 1;
end

% enter debug mode here
if I.keyboard
    keyboard;
end

% string with modeling parameters
always_include = {'distr'};
always_exclude = {...
    'run', 'figh', 'keyboard', 'plot_figure', 'plot_win', ...
    'plotcausal', 'overwrite'};
param_string_modelfit = optInputs_to_string(I, C_value, always_include, always_exclude);

% file to save results to
% can optionally and just return this file
MAT_file = mkpdir([L.output_directory '/model_fit_' param_string_modelfit '.mat']);
if ~I.run
    M = struct;
    return;
end

% function used to transform predictions
switch I.tranpred
    case 'none'
        tranpred = @(x)x;
    case 'sqrt'
        tranpred = @(x)sqrt(x);
    case 'sq'
        tranpred = @(x)x.^2;
    otherwise
        error('No matching tranpred parameter');
end

% function used to transform correlation
switch I.tranweightnsegs
    case 'none'
        tranweightnsegsfn = @(x)x;
    case 'sqrt'
        tranweightnsegsfn = @(x)sqrt(x);
    case 'sqrtm3'
        tranweightnsegsfn = @(x)sqrt(x-3);
    otherwise
        error('No matching tranweightnsegs parameter');
end

% function used to transform weights
switch I.tranweightdenom
    case 'none'
        tranweightdenomfn = @(x)x;
    case 'square'
        tranweightdenomfn = @(x)x.^2;
    otherwise
        error('No matching tranweightdenom parameter');
end

%% Model fit

if ~exist(MAT_file, 'file') || I.overwrite
    
    [n_lags, n_seg, n_channels] = size(L.same_context);
    assert(n_channels == length(L.channels));
    assert(length(L.unique_segs)==n_seg);
    assert(n_lags == length(L.lag_t))
    
    % measure to predict with the model
    same_context = L.same_context;
    diff_context = L.diff_context;
    denom_factor = nan(size(L.same_context));
    M.Y = nan(size(same_context));
    for k = 1:n_seg
        if I.normcorr
            assert(~I.divnorm);
            baseline = mean(same_context(:,k,:));
            M.Y(:,k,:) = 1 - bsxfun(@times, (same_context(:,k,:)-diff_context(:,k,:)), 1./baseline);
            denom_factor(:,k,:) = repmat(baseline, 1, size(same_context,1));
        elseif I.divnorm
            assert(~I.normcorr);
            for q = 1:size(diff_context,3)
                xi = same_context(:,k,q)>I.mindenom;
                M.Y(xi,k,q) = diff_context(xi,k,q) ./ same_context(xi,k,q);
                denom_factor(xi,k,q) = same_context(xi,k,q);
                clear xi;
            end
        else
            M.Y(:,k,:) = same_context(:,k,:)-diff_context(:,k,:);
            denom_factor(:,k,:) = 1;
        end
    end
    
    % segment dependent weights
    w_segs = tranweightnsegsfn(L.n_total_segs);
    
    % valid segment durations
    valid_seg_durs = find(w_segs>0);
        
    % integration period
    M.intper_sec = logspace(log10(I.intper_range(1)), log10(I.intper_range(2)), I.nintper);
    M.delay_sec = I.delay_range(1):I.delay_interval:I.delay_range(2);
    M.shape = I.shape;
    M.err = nan(length(M.intper_sec), length(M.delay_sec), length(M.shape), n_channels);
    % M.Yh = nan([n_lags, n_seg, length(M.intper_sec), length(M.delay_sec), length(M.shape), n_channels]);
    M.causal_win = false(length(M.intper_sec), length(M.delay_sec), length(M.shape));
    for m = 1:length(M.shape)
        for i = 1:length(M.intper_sec)
            fprintf('shape %.2f, intper %.0f ms\n', M.shape(m), M.intper_sec(i)*1000);
            drawnow;
            for j = 1:length(M.delay_sec)
                
                Yh = nan(n_lags, n_seg, n_channels);
                
                % calculate predictions for each segment
                for k = valid_seg_durs
                    
                    [winpow, ~, overlap] = win_power_ratio(L.unique_segs(k)/1000, ...
                        I.distr, M.intper_sec(i), M.delay_sec(j), ...
                        'shape', M.shape(m), 'forcecausal', I.forcecausal, ...
                        'tsec', L.lag_t, 'rampwin', I.rampwin, 'rampdur', I.rampdur);
                    if I.winpowratio
                        predictor_notrans = winpow;
                    else
                        predictor_notrans = overlap;
                    end
                    predictor_notrans(predictor_notrans<0) = 0;
                    predictor = tranpred(predictor_notrans);
                    
                    % prdictions for each electrode
                    for q = 1:n_channels
                        if I.normcorr || I.divnorm
                            Yh(:,k,q) = predictor;
                        else
                            Yh(:,k,q) = (1-predictor)*pinv(1-predictor)*M.Y(:,k,q);
                        end
                    end
                end
                
                % calculate error
                for q = 1:n_channels
                    E = (Yh(:,valid_seg_durs,q) - M.Y(:,valid_seg_durs,q)).^2;
                    if I.weightdenom
                        w_denom = denom_factor(:,valid_seg_durs,q);
                        w_denom(~(w_denom>I.mindenom)) = 0;
                        w_denom = tranweightdenomfn(w_denom);
                        w_total = bsxfun(@times, w_denom, w_segs(valid_seg_durs));
                    else
                        w_total = w_segs(valid_seg_durs);
                    end
                    w_total = w_total / sum(w_total(:));
                    Ew = bsxfun(@times, E, w_total);
                    clear w_total w_denom;
                    M.err(i,j,m,q) = nanmean(Ew(:));
                end
            end
        end
    end
    
    % find best prediction
    M.Ybest = nan(n_lags, n_seg, n_channels);
    M.best_intper_sec = nan(n_channels,1);
    M.best_delay_sec = nan(n_channels,1);
    M.best_shape_smp = nan(n_channels,1);
    for q = 1:n_channels
        X = M.err(:,:,:,q);
        if any(M.causal_win(:)) && I.bestcausal
            X(~M.causal_win) = inf;
        end
        [~,xi] = min(X(:));
        [a,b,c] = ind2sub(size(X), xi);
        % M.Ybest(:,:,q) = M.Yh(:,:,a,b,c,q);
        M.best_intper_sec(q) = M.intper_sec(a);
        M.best_delay_sec(q) = M.delay_sec(b);
        M.best_shape(q) = M.shape(c);
        clear X;
        
        % get prediction
        for k = valid_seg_durs
            
            [winpow, ~, overlap] = win_power_ratio(L.unique_segs(k)/1000, I.distr, ...
                M.best_intper_sec(q), M.best_delay_sec(q), ...
                'shape', M.best_shape(q), 'forcecausal', I.forcecausal, ...
                'tsec', L.lag_t, 'rampwin', I.rampwin, 'rampdur', I.rampdur);
            if I.winpowratio
                predictor_notrans = winpow;
            else
                predictor_notrans = overlap;
            end
            predictor_notrans(predictor_notrans<0) = 0;
            predictor = tranpred(predictor_notrans);
            
            if I.normcorr || I.divnorm
                M.Ybest(:,k,q) = predictor;
            else
                M.Ybest(:,k,q) = (1-predictor)*pinv(1-predictor)*M.Y(:,k,q);
            end
        end
    end
    
    M.channels = L.channels;
    
    save(MAT_file, 'M', '-v7.3');
    
else
    
    load(MAT_file, 'M')
    
end

%% Plot the results

if I.plot_figure
        
    % create a figure handle if not supplied
    if isempty(I.figh)
        I.figh = figure;
    end
    
    plot_win_string = [num2str(I.plot_win(1)) '-' num2str(I.plot_win(2))];
    
    n_channels = length(L.channels);
    for q = 1:n_channels
        chan = L.channels(q);
        
        if strcmp(L.chnames{chan}, ['ch' num2str(chan)])
            chname = L.chnames{chan};
        else
            chname = ['ch' num2str(chan) '-' L.chnames{chan}];
        end
        
        % plot prediction for best delay, lineplot
        clf(I.figh);
        set(I.figh, 'Position', [100, 100, 900, 900]);
        if I.normcorr || I.divnorm
            corr_range = [-0.5 1.5];
            invariance_line = 1;
        else
            X = M.Y(:,:,q);
            corr_range = quantile(X(:), [0.01, 0.99]);
            clear X;
            invariance_line = NaN;
        end
        valid_seg_durs = find(L.n_total_segs>0);
        for k = valid_seg_durs
            subplot(4, 2, k);
            if I.normcorr
                plot(L.lag_t([1 end]) * 1000, [1,1], 'k--', 'LineWidth', 2);
            end
            hold on;
            plot(L.lag_t([1 end]) * 1000, [0,0], 'k--', 'LineWidth', 2);
            h1 = plot(L.lag_t * 1000, M.Y(:,k,q), 'LineWidth', 2);
            h2 = plot(L.lag_t * 1000, M.Ybest(:,k,q), 'LineWidth', 2);
            plot(L.unique_segs(k)*[1 1], corr_range, 'k--', 'LineWidth', 2);
            if ~isnan(invariance_line); plot(I.plot_win * 1000, invariance_line*[1 1], 'k--', 'LineWidth', 2); end
            xlim(I.plot_win * 1000);
            ylim(corr_range);
            xlabel('Time Lag (ms)');
            ylabel('Pearson Correlation');
            title(sprintf('Seg: %.0f ms', L.unique_segs(k)))
            if k == 1
                legend([h1, h2], 'Data', 'Model');
            end
            box off;
        end
        fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-prediction-lineplot_plotwin' plot_win_string]);
        print_wrapper([fname '.png']);
        print_wrapper([fname '.pdf']);
        savefig(I.figh, mkpdir([fname '.fig']));
        
        % plot prediction for best delay, image
        clf(I.figh);
        set(I.figh, 'Position', [100, 100, 900, 600]);
        subplot(2,1,1);
        imagesc(M.Y(:,:,q)', corr_range(2) * [-1, 1]);
        subplot(2,1,2);
        imagesc(M.Ybest(:,:,q)', corr_range(2) * [-1, 1]);
        for i = 1:2
            subplot(2,1,i);
            colorbar;
            colormap(flipud(cbrewer('div', 'RdBu', 128)));
            set(gca, 'YTick', 1:length(L.unique_segs), 'YTickLabel', L.unique_segs);
            xticks = get(gca, 'XTick');
            set(gca, 'XTick', xticks, 'XTickLabel', L.lag_t(xticks)*1000);
            ylabel('Seg Dur (ms)'); xlabel('Lag (ms)');
            if i == 2
                title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q)*1000, M.best_delay_sec(q)*1000));
            end
        end
        fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-prediction_plotwin' plot_win_string]);
        print_wrapper([fname '.png']);
        savefig(I.figh, mkpdir([fname '.fig']));
        % export_fig([fname '.png'], '-png', '-transparent', '-r100');
        
        % plot the error vs. parameters
        X = M.err(:,:,M.best_shape(q)==M.shape,q);
        if any(M.causal_win(:)) && I.plotcausal
            causal_string = '_only-causal';
            X(~M.causal_win(:,:,M.best_shape(q)==M.shape)) = NaN;
        else
            causal_string = '';
        end
        bounds = quantile(-X(:), [0.8, 1]);
        if ~all(isnan(X(:)))
            clf(I.figh);
            set(I.figh, 'Position', [100, 100, 600, 600]);
            
            imagesc(-X, bounds);
            colormap(parula(128));
            colorbar;
            yticks = round(get(gca, 'YTick'));
            set(gca, 'YTick', yticks, 'YTickLabel', 1000*M.intper_sec(yticks));
            % ytick = interp1(log2(M.intper_sec), 1:length(M.intper_sec), log2(M.intper_sec));
            % set(gca, 'YTick', ytick, 'YTickLabel', M.intper_sec*1000);
            xticks = round(get(gca, 'XTick'));
            set(gca, 'XTick', xticks, 'XTickLabel', 1000*M.delay_sec(xticks));
            xlabel('Delay (ms)'); ylabel('Receptive Field Duration (ms)');
            title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q)*1000, M.best_delay_sec(q)*1000));
            set(gca, 'FontSize', 12);
            fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-error' causal_string]);
            print_wrapper([fname '.png']);
            savefig(I.figh, mkpdir([fname '.fig']));
            % export_fig([fname '.png'], '-png', '-transparent', '-r100');
            clear X;
        end
    end
    
end