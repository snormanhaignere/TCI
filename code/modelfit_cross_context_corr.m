function [M, MAT_file] = modelfit_cross_context_corr(L, varargin)

% Fit parametric window using cross-context correlation data. This function
% should be applied to the output structure L returned by
% cross_context_corr.m

%% Default parameters

clear I;

% loss function, options are:
% sqerr: squared error
% nse: normalized squared error
% corr: negative Pearson correlation
I.lossfn = 'sqerr';

% distribution used to model the window
% 'gamma' or 'gauss'
I.distr = 'gamma';

% shape parameter for parametric window, not relevant for gaussian
% 1 -> exponential, higher -> more gaussian
I.shape = 1;

% central interval used to calculate integration period
I.centralinterval = 0.75;

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
I.shape = [1,2,3,5,10]; 

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

% whether to create null samples with phase scrambling
I.nullsmps = 0;

% seed to create fixed random decision (e.g. for phase scrambling)
I.randseed = 1;

% whether to overwrite existing results
I.overwrite = false;

% whether to enter debug mode
I.keyboard = false;

% whether to plot figures
I.plot_figure = true;

% window used for plotting
I.plot_win = L.lag_t([1 end]); % in seconds

% quantile of error map to plot
I.ploterrquant = 0.1;

% whether to plot bootstrapped samples
I.plot_extrasmps = false;

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
    'plotcausal', 'overwrite', 'ploterrquant', 'plot_extrasmps'};
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

% function used to evaluate lossfn
switch I.lossfn
    case 'sqerr'
        lossfn = @(x,y,w)weighted_squared_error(x,y,w);
    case 'nse'
        lossfn = @(x,y,w)weighted_nse(x,y,w);
    case 'corr'
        lossfn = @(x,y,w)(-weighted_pearson_corr(x,y,w));
    otherwise
        error('No matching loss for %s', I.lossfn);
end

[n_lags, n_seg, n_channels, ~] = size(L.same_context);
assert(n_channels == length(L.channels));
assert(length(L.unique_segs)==n_seg);
assert(n_lags == length(L.lag_t))

%% Model fit

if ~exist(MAT_file, 'file') || I.overwrite
    
    ResetRandStream2(I.randseed);
    
    % can optionally create null samples via phase scrambling
    if I.nullsmps==0
        same_context = L.same_context;
        diff_context = L.diff_context;
        denom_factor = nan(size(L.same_context));
    else
        assert(size(L.diff_context,4)==1); % there should only be one sample already (i.e. no bootstrapping)
        diff_context_null = nan([n_lags, n_seg, n_channels, I.nullsmps]);
        for i = 1:I.nullsmps
            diff_context_null(:,:,:,i) = phasescram(L.diff_context);
        end
        diff_context = cat(4, L.diff_context, diff_context_null);
        same_context = repmat(L.same_context, [1, 1, 1, I.nullsmps+1]);
        denom_factor = nan(size(diff_context));
    end
    n_smps = size(same_context,4);

    % create the vector to be predicted
    M.Y = nan(size(same_context));
    if I.normcorr
        assert(~I.divnorm);
        baseline = mean(same_context);
        M.Y = 1 - bsxfun(@times, (same_context-diff_context), 1./baseline);
        denom_factor = repmat(baseline, [size(same_context,1), ones(1, ndims(same_context)-1)]);
    elseif I.divnorm
        assert(~I.normcorr);
        xi = same_context > I.mindenom;
        M.Y(xi) = diff_context(xi) ./ same_context(xi);
        denom_factor(xi) = same_context(xi);
        clear xi;
    else
        M.Y = same_context-diff_context;
        denom_factor(:) = 1;
    end
        
    % valid segment durations
    valid_seg_durs = find(L.n_total_segs(:,1)>0)';
    
    % integration period
    M.intper_sec = logspace(log10(I.intper_range(1)), log10(I.intper_range(2)), I.nintper);
    M.delay_sec = I.delay_range(1):I.delay_interval:I.delay_range(2);
    M.shape = I.shape;
    M.loss = nan(length(M.intper_sec), length(M.delay_sec), length(M.shape), n_channels, n_smps);
    % M.Yh = nan([n_lags, n_seg, length(M.intper_sec), length(M.delay_sec), length(M.shape), n_channels]);
    M.causal_win = false(length(M.intper_sec), length(M.delay_sec), length(M.shape));
    for m = 1:length(M.shape)
        for i = 1:length(M.intper_sec)
            fprintf('shape %.2f, intper %.0f ms\n', M.shape(m), M.intper_sec(i)*1000);
            drawnow;
            for j = 1:length(M.delay_sec)
                
                Y_model = nan(n_lags, n_seg);
                
                % calculate predictions for each segment
                for k = valid_seg_durs
                    
                    [winpow, ~, overlap, M.causal_win(i,j,m)] = win_power_ratio(L.unique_segs(k)/1000, ...
                        I.distr, M.intper_sec(i), M.delay_sec(j), ...
                        'shape', M.shape(m), 'forcecausal', I.forcecausal, ...
                        'tsec', L.lag_t, 'rampwin', I.rampwin, 'rampdur', I.rampdur, ...
                        'centralinterval', I.centralinterval, 'delaypoint', I.delaypoint);
                    if I.winpowratio
                        predictor_notrans = winpow;
                    else
                        predictor_notrans = overlap;
                    end
                    predictor_notrans(predictor_notrans<0) = 0;
                    predictor = tranpred(predictor_notrans);
                    
                    Y_model(:,k) = predictor;
                    
                end
                
                if I.normcorr || I.divnorm
                    Y_pred = repmat(Y_model, [1, 1, n_channels, n_smps]);
                else
                    Y_pred = nan(size(same_context));
                    for s = 1:n_smps
                        for q = 1:n_channels
                            for k = valid_seg_durs
                                Y_pred(:,k,q,s) = (1-predictor)*pinv(1-predictor)*M.Y(:,k,q,s);
                            end
                        end
                    end
                end
                
                % denominator weights
                % lag x segdur x channel x sample
                W_denom = denom_factor;
                W_denom(~(W_denom>I.mindenom)) = 0;
                W_denom = tranweightdenomfn(W_denom);
                
                % segment dependent weights
                % segdur x sample (i.e. for bootstrapping)
                W_segs = tranweightnsegsfn(L.n_total_segs);
                
                % combine
                W_total = bsxfun(@times, W_denom, reshape(W_segs, [1, size(W_segs,1), 1, size(W_segs,2)]));
                
                % format
                Y_pred_valid = Y_pred(:,valid_seg_durs,:,:);
                Y_data_valid = M.Y(:,valid_seg_durs,:,:);
                W_total_valid = W_total(:,valid_seg_durs,:,:);
                Y_pred_valid = reshape(Y_pred_valid, [n_lags * length(valid_seg_durs), n_channels, n_smps]);
                Y_data_valid = reshape(Y_data_valid, [n_lags * length(valid_seg_durs), n_channels, n_smps]);
                W_total_valid = reshape(W_total_valid, [n_lags * length(valid_seg_durs), n_channels, n_smps]);
                
                % calculate loss
                M.loss(i,j,m,:,:) = lossfn(Y_pred_valid, Y_data_valid, W_total_valid);
                
                %                 % calculate error
                %                 for q = 1:n_channels
                %                     %E = (Yh(:,valid_seg_durs,q) - M.Y(:,valid_seg_durs,q)).^2;
                %                     if I.weightdenom
                %                         W_denom = denom_factor(:,valid_seg_durs,q);
                %                         W_denom(~(W_denom>I.mindenom)) = 0;
                %                         W_denom = tranweightdenomfn(W_denom);
                %                     else
                %                         W_total = W_segs(valid_seg_durs);
                %                     end
                %                     % w_total = w_total / sum(w_total(:));
                %                     X1 = Y_pred(:,valid_seg_durs,q);
                %                     X2 = M.Y(:,valid_seg_durs,q);
                %                     % Ew = bsxfun(@times, E, w_total);
                %                     M.loss(i,j,m,q) = lossfn(X1(:), X2(:), W_total(:));
                %                     clear w_total w_denom;
                %                 end
            end
        end
    end
    
    % find best prediction
    M.Ybest = nan(n_lags, n_seg, n_channels, n_smps);
    M.best_intper_sec = nan(n_channels,n_smps);
    M.best_delay_sec = nan(n_channels,n_smps);
    M.best_shape = nan(n_channels,n_smps);
    M.best_loss = nan(n_channels,n_smps);
    for s = 1:n_smps
        for q = 1:n_channels
            X = M.loss(:,:,:,q,s);
            if any(M.causal_win(:)) && I.bestcausal
                X(~M.causal_win) = inf;
            end
            [M.best_loss(q,s),xi] = min(X(:));
            [a,b,c] = ind2sub(size(X), xi);
            M.best_intper_sec(q,s) = M.intper_sec(a);
            M.best_delay_sec(q,s) = M.delay_sec(b);
            M.best_shape(q,s) = M.shape(c);
            clear X;
            
            % get prediction
            for k = valid_seg_durs
                
                [winpow, ~, overlap] = win_power_ratio(L.unique_segs(k)/1000, I.distr, ...
                    M.best_intper_sec(q,s), M.best_delay_sec(q,s), ...
                    'shape', M.best_shape(q,s), 'forcecausal', I.forcecausal, ...
                    'tsec', L.lag_t, 'rampwin', I.rampwin, 'rampdur', I.rampdur, ...
                    'centralinterval', I.centralinterval, 'delaypoint', I.delaypoint);
                if I.winpowratio
                    predictor_notrans = winpow;
                else
                    predictor_notrans = overlap;
                end
                predictor_notrans(predictor_notrans<0) = 0;
                predictor = tranpred(predictor_notrans);
                
                if I.normcorr || I.divnorm
                    M.Ybest(:,k,q,s) = predictor;
                else
                    M.Ybest(:,k,q,s) = (1-predictor)*pinv(1-predictor)*M.Y(:,k,q,s);
                end
            end
        end
    end
    
    if I.nullsmps>0
        M.logP_gaussfit = sig_via_null_gaussfit(-M.best_loss(:,1)', -M.best_loss(:,2:end)');
        M.logP_counts = sig_via_null_counts(-M.best_loss(:,1)', -M.best_loss(:,2:end)');
    end
    
    M.channels = L.channels;
    
    save(MAT_file, 'M', '-v7.3');
    
else
    
    load(MAT_file, 'M')
    
end

M.param_string_modelfit = param_string_modelfit;

%% Plot the results

if I.plot_figure
    
    % create a figure handle if not supplied
    if isempty(I.figh)
        I.figh = figure;
    end
    
    plot_win_string = [num2str(I.plot_win(1)) '-' num2str(I.plot_win(2))];
    
    n_smps = size(M.Y,4);
    n_channels = length(L.channels);
    if I.plot_extrasmps
        smps = 1:n_smps;
    else
        smps = 1;
    end
    for q = 1:n_channels
        for s = smps
            chan = L.channels(q);
            
            if strcmp(L.chnames{chan}, ['ch' num2str(chan)])
                chname = L.chnames{chan};
            else
                chname = ['ch' num2str(chan) '-' L.chnames{chan}];
            end
            
            if s > 1
                smpstr = ['_smp' num2str(s-1)];
            else
                smpstr = '';
            end
            
            % plot prediction for best delay, lineplot
            clf(I.figh);
            set(I.figh, 'Position', [100, 100, 900, 900]);
            if I.normcorr || I.divnorm
                corr_range = [-0.5 1.5];
                invariance_line = 1;
            else
                X = M.Y(:,:,q,s);
                corr_range = quantile(X(:), [0.01, 0.99]);
                clear X;
                invariance_line = NaN;
            end
            valid_seg_durs = find(L.n_total_segs(:,1)>0)';
            for k = valid_seg_durs
                subplot(4, 2, k);
                if I.normcorr
                    plot(L.lag_t([1 end]) * 1000, [1,1], 'k--', 'LineWidth', 2);
                end
                hold on;
                plot(L.lag_t([1 end]) * 1000, [0,0], 'k--', 'LineWidth', 2);
                h1 = plot(L.lag_t * 1000, M.Y(:,k,q,s), 'LineWidth', 2);
                h2 = plot(L.lag_t * 1000, M.Ybest(:,k,q,s), 'LineWidth', 2);
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
            fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-prediction-lineplot_plotwin' plot_win_string smpstr]);
            print_wrapper([fname '.png']);
            print_wrapper([fname '.pdf']);
            savefig(I.figh, mkpdir([fname '.fig']));
            
            % plot prediction for best delay, image
            clf(I.figh);
            set(I.figh, 'Position', [100, 100, 900, 600]);
            subplot(2,1,1);
            imagesc(M.Y(:,:,q,s)', corr_range(2) * [-1, 1]);
            subplot(2,1,2);
            imagesc(M.Ybest(:,:,q,s)', corr_range(2) * [-1, 1]);
            for i = 1:2
                subplot(2,1,i);
                colorbar;
                colormap(flipud(cbrewer('div', 'RdBu', 128)));
                set(gca, 'YTick', 1:length(L.unique_segs), 'YTickLabel', L.unique_segs);
                xticks = get(gca, 'XTick');
                set(gca, 'XTick', xticks, 'XTickLabel', L.lag_t(xticks)*1000);
                ylabel('Seg Dur (ms)'); xlabel('Lag (ms)');
                if i == 2
                    title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q,s)*1000, M.best_delay_sec(q,s)*1000));
                end
            end
            fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-prediction_plotwin' plot_win_string smpstr]);
            print_wrapper([fname '.png']);
            savefig(I.figh, mkpdir([fname '.fig']));
            % export_fig([fname '.png'], '-png', '-transparent', '-r100');
            
            % plot the error vs. parameters
            X = M.loss(:,:,M.best_shape(q,s)==M.shape,q,s);
            if any(M.causal_win(:)) && I.plotcausal
                causal_string = '_only-causal';
                X(~M.causal_win(:,:,M.best_shape(q,s)==M.shape)) = inf;
            else
                causal_string = '';
            end
            
            cmap = flipud(cbrewer('seq', 'Reds', 128));
            [minX,zi] = min(X(:));
            [xi, ~] = ind2sub(size(X), zi);
            bounds = [minX, quantile(X(xi,:), I.ploterrquant)];
            clear xi zi;
            % bounds = quantile(-X(:), [1-I.ploterrquant, 1]);
            if ~all(isnan(X(:)))
                clf(I.figh);
                set(I.figh, 'Position', [100, 100, 600, 600]);
                
                imagesc(X, bounds);
                colormap(cmap);
                colorbar;
                yticks = round(get(gca, 'YTick'));
                set(gca, 'YTick', yticks, 'YTickLabel', 1000*M.intper_sec(yticks));
                % ytick = interp1(log2(M.intper_sec), 1:length(M.intper_sec), log2(M.intper_sec));
                % set(gca, 'YTick', ytick, 'YTickLabel', M.intper_sec*1000);
                xticks = round(get(gca, 'XTick'));
                set(gca, 'XTick', xticks, 'XTickLabel', 1000*M.delay_sec(xticks));
                xlabel('Delay (ms)'); ylabel('Receptive Field Duration (ms)');
                title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q,s)*1000, M.best_delay_sec(q,s)*1000));
                set(gca, 'FontSize', 12);
                fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-error' causal_string smpstr]);
                print_wrapper([fname '.png']);
                savefig(I.figh, mkpdir([fname '.fig']));
                % export_fig([fname '.png'], '-png', '-transparent', '-r100');
                % clear X;
            end
        end
    end
    
end