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
I.lossfn = 'unbiased-sqerr';

% how to estimate error if using unbiased squared error
I.esterr = 'splithalf-mean';

% distribution used to model the window
% only 'gamma' allowed currently
I.distr = 'gamma';

% shape parameter for parametric window, not relevant for gaussian
% 1 -> exponential, higher -> more gaussian
I.shape = 1;

% density interval used to calculate integration period
% intervaltype:
% 'highest': highest density interval (default, the minimal interval of a given mass)
% 'center': central interval, i.e. from [0+(1-intervalmass)/2, 1-(1-intervalmass)/2]
% 'start': starting interval, i.e. from [0 to intervalmass]
I.intervaltype = 'highest';
I.intervalmass = 0.75;

% range of integration periods to consider in seconds
I.intper_range = L.unique_segs([1 end])/1000; % in seconds

% number of integration periods to consider between this range
% spaced on a logarithmic scale
I.nintper = 20;

% range of delays to consider in seconds
I.delay_range = [0, L.unique_segs(end)/1000];

% interval between delays in seconds
I.delay_interval = 1/L.sr;

% whether and by how much to ramp segments
% this should match the actual ramps used
% in the stitching together of segments
I.rampwin = L.rampwin;
I.rampdur = L.rampdur; % in seconds

% lags with a same-context correlation below this value are excluded
I.minreliability = 0.01;

% range of shape parameters for gamma distribution
% 1 -> exponential, higher -> more Gaussian
I.shape = [1,2,3,5];

% whether to use the power ratio between
% segments to predict the correlation
I.winpowratio = true;

% whether to transform the weighting applied to segments
I.tranweightnsegs = 'none'; % applied to total segs

% whether to transform the predictions
I.tranpred = 'none';

% whether to create null samples with phase scrambling
I.nullsmps = 0;

% skip calculation of the best predictors for null samples
I.skipnullpreds = true;

% skip calculation of best predictors for extra samples
I.skipsmppreds = ~isfield(L, 'splitsmps');

% seed to create fixed random decision (e.g. for phase scrambling)
I.randseed = 1;

% whether to overwrite existing results
I.overwrite = false;

% whether to enter debug mode
I.keyboard = false;

% how to quantify delay in the plots
I.plot_delaystat = 'median';

% range of delays to plot
I.plot_delay_range = [0, L.unique_segs(end)/1000];

% whether to plot figures
I.plot_figure = true;

% window used for plotting
I.plot_win = L.lag_t([1 end]); % in seconds

% line width for plotting
I.linewidth = 2;

% quantile of error map to plot
I.ploterrquant = 0.3;

% whether to plot bootstrapped samples
I.plot_extrasmps = isfield(L, 'splitsmps');
I.plot_nullsmps = false;

% whether to run the analysis
% if false, just returns the parameters to be used
I.run = true;

% figure handle used for plotting
I.figh = matlab.ui.Figure.empty;

%% Parse user-specified parameters, create parameter string for this analysis

[I, C, C_value] = parse_optInputs_keyvalue(varargin, I);

% enter debug mode here
if I.keyboard
    keyboard;
end

% set plotting range to actual range, if not otherwise specified
if C.delay_range && ~C.plot_delay_range
    I.plot_delay_range = I.delay_range;
end

% string with modeling parameters
always_include = {'lossfn'};
always_exclude = {...
    'run', 'figh', 'keyboard', 'plot_figure', 'plot_win', 'linewidth', ...
    'overwrite', 'ploterrquant', 'plot_extrasmps', 'plot_delaystat', 'plot_delay_range'};
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

% function used to evaluate lossfn
switch I.lossfn
    case 'sqerr'
        lossfn = @(x,y,w)sum(w.*bsxfun(@minus, x, y).^2,1);
    case 'unbiased-sqerr'
        lossfn = @(x,y,w,e)sum(w.*bsxfun(@minus, x, y).^2 - w.*e,1);
    case 'nse'
        lossfn = @(x,y,w)weighted_nse(x,y,w,true);
    case 'nse-subtract'
        lossfn = @(x,y,w)weighted_nse(x,y,w,false);
    case 'corr'
        lossfn = @(x,y,w)(-weighted_pearson_corr(x,y,w));
    otherwise
        error('No matching loss for %s', I.lossfn);
end

[n_lags, n_seg, n_channels, n_smps] = size(L.same_context);
assert(n_channels == length(L.channels));
assert(length(L.unique_segs)==n_seg);
assert(n_lags == length(L.lag_t))

%% Model fit

if ~exist(MAT_file, 'file') || I.overwrite
    
    clear M;
    ResetRandStream2(I.randseed);
    
    % error to use for noise correction
    if strcmp(I.lossfn, 'unbiased-sqerr')
        switch I.esterr
            case 'bstrap'
                M.same_context_err = L.same_context_bstrap_err;
            case 'splithalf'
                M.same_context_err = L.same_context_err;
            case 'splithalf-mean'
                M.same_context_err = mean(L.same_context_err,1);
                M.same_context_err = repmat(M.same_context_err, [n_lags, ones(1,ndims(L.same_context_err)-1)]);
            case 'neglag'
                M.same_context_err = mean(L.diff_context(L.lag_t<0,:,:,:).^2,1);
                M.same_context_err = repmat(M.same_context_err, [n_lags, ones(1,ndims(L.diff_context)-1)]);
            case 'neglag-smps'
                xi = randi(sum(L.lag_t<0), [n_lags,1]);
                M.same_context_err = L.diff_context(xi, :, :, :).^2;
            otherwise
                error('No matching case for esterr=%s\n', I.esterr);
        end
    end
    
    % can optionally create null samples via phase scrambling
    if I.nullsmps==0
        M.same_context = L.same_context;
        M.diff_context = L.diff_context;
    else
        diff_context_null = nan([n_lags, n_seg, n_channels, n_smps, I.nullsmps]);
        for i = 1:I.nullsmps
            diff_context_null(:,:,:,:,i) = real(phasescram(L.diff_context));
        end
        M.diff_context = cat(5, L.diff_context, diff_context_null);
        M.same_context = repmat(L.same_context, [1, 1, 1, 1, I.nullsmps+1]);
        M.same_context_err = repmat(M.same_context_err, [1, 1, 1, 1, I.nullsmps+1]);
    end
    
    % valid segment durations
    valid_seg_durs = find(L.n_total_segs(:,1)>0)';
    n_valid_segdurs = length(valid_seg_durs);
    
    % segment dependent weights
    % segdur x sample (i.e. for bootstrapping)
    W_segs = tranweightnsegsfn(L.n_total_segs);
    
    % reshape
    % lag x segment duration x channel x sample x null sample
    W_total = reshape(W_segs, [1, size(W_segs,1), 1, size(W_segs,2)]);
    W_total = repmat(W_total, [n_lags, 1, n_channels, 1, I.nullsmps+1]);
    
    % select valid segments
    same_context_valid = M.same_context(:,valid_seg_durs,:,:,:);
    diff_context_valid = M.diff_context(:,valid_seg_durs,:,:,:);
    W_total_valid = W_total(:,valid_seg_durs,:,:,:);
    if strcmp(I.lossfn, 'unbiased-sqerr')
        same_context_err_valid = M.same_context_err(:,valid_seg_durs,:,:,:);
    end
    
    % unwrap lag and segment duration into one vectory
    same_context_format = reshape(same_context_valid, [n_lags * length(valid_seg_durs), n_channels, n_smps, I.nullsmps+1]);
    diff_context_format = reshape(diff_context_valid, [n_lags * length(valid_seg_durs), n_channels, n_smps, I.nullsmps+1]);
    W_total_format = reshape(W_total_valid, [n_lags * length(valid_seg_durs), n_channels, n_smps, I.nullsmps+1]);
    if strcmp(I.lossfn, 'unbiased-sqerr')
        same_context_err_format = reshape(same_context_err_valid, [n_lags * length(valid_seg_durs), n_channels, n_smps, I.nullsmps+1]);
    end
    
    % set weights for invalid lags to zero
    invalid_lags = same_context_format < I.minreliability;
    W_total_format(invalid_lags) = 0;
    same_context_format(invalid_lags) = 0;
    diff_context_format(invalid_lags) = 0;
    if strcmp(I.lossfn, 'unbiased-sqerr')
        same_context_err_format(invalid_lags) = 0;
    end
    
    % normalize weights
    W_total_format = bsxfun(@times, W_total_format, 1./sum(W_total_format,1));
    
    % integration period
    M.intper_sec = logspace(log10(I.intper_range(1)), log10(I.intper_range(2)), I.nintper);
    M.delay_sec_start = I.delay_range(1):I.delay_interval:I.delay_range(2);
    M.delay_sec_start = round(M.delay_sec_start*L.sr)/L.sr;
    M.shape = I.shape;
    M.loss = nan(length(M.intper_sec), length(M.delay_sec_start), length(M.shape), n_channels, n_smps, I.nullsmps+1);
    for m = 1:length(M.shape)
        for i = 1:length(M.intper_sec)
            fprintf('shape %.2f, intper %.0f ms\n', M.shape(m), M.intper_sec(i)*1000);drawnow;
            
            % calculate predictions for each segment
            Y_model = nan(n_lags, n_valid_segdurs);
            for k = 1:n_valid_segdurs
                [winpow, ~, overlap] = win_power_ratio(L.unique_segs(valid_seg_durs(k))/1000, ...
                    I.distr, M.intper_sec(i), 0, ...
                    'shape', M.shape(m), 'tsec', L.lag_t, ...
                    'rampwin', I.rampwin, 'rampdur', I.rampdur, ...
                    'intervalmass', I.intervalmass, ...
                    'intervaltype', I.intervaltype, ...
                    'delaypoint', 'start');
                if I.winpowratio
                    predictor_notrans = winpow;
                else
                    predictor_notrans = overlap;
                end
                predictor_notrans(predictor_notrans<0) = 0;
                predictor = tranpred(predictor_notrans);
                
                Y_model(:,k) = predictor;
            end
            
            % create delays
            % lag x segment duration x time delays
            Y_model_with_delays = add_delays(Y_model, checkint(M.delay_sec_start*L.sr));
            
            for j = 1:length(M.delay_sec_start)
                
                % multiply model prediction by same-context correlation
                % to get a prediction of the different context correlation
                X = Y_model_with_delays(:,:,j);
                diff_context_pred = bsxfun(@times, same_context_format, X(:));
                
                % calculate loss
                if strcmp(I.lossfn, 'unbiased-sqerr')
                    err = bsxfun(@times, same_context_err_format, X(:).^2);
                    M.loss(i,j,m,:,:,:) = lossfn(diff_context_format, diff_context_pred, W_total_format, err);
                else
                    M.loss(i,j,m,:,:,:) = lossfn(diff_context_format, diff_context_pred, W_total_format);
                end
            end
        end
    end
    
    % find best prediction
    M.diff_context_bestpred = nan(n_lags, n_seg, n_channels, (n_smps-1)*I.skipsmppreds+1, I.nullsmps*I.skipnullpreds+1);
    M.best_intper_sec = nan(n_channels, n_smps, I.nullsmps+1);
    M.best_delay_sec_start = nan(n_channels, n_smps, I.nullsmps+1);
    M.best_delay_sec_median = nan(n_channels, n_smps, I.nullsmps+1);
    M.best_shape = nan(n_channels, n_smps, I.nullsmps+1);
    M.best_loss = nan(n_channels, n_smps, I.nullsmps+1);
    for q = 1:n_channels
        for s = 1:n_smps
            
            % if using cross validation
            % compute the best model using the mean of
            % all other samples (from non-permuted data)
            if isfield(L, 'splitsmps') && any(s == L.splitsmps) 
                X = mean(M.loss(:,:,:,q,setdiff(L.splitsmps, s),1),4);
                [~,xi] = min(X(:));
                [a,b,c] = ind2sub(size(X), xi);
                clear X;
                splitsmp = true;
            else
                splitsmp = false;
            end
            
            for l = 1:I.nullsmps+1
                if ~splitsmp
                    X = M.loss(:,:,:,q,s,l);
                    [~,xi] = min(X(:));
                    [a,b,c] = ind2sub(size(X), xi);
                    splitsmp = true;
                    clear X;
                end
                M.best_loss(q,s,l) = M.loss(a,b,c,q,s,l);
                M.best_intper_sec(q,s,l) = M.intper_sec(a);
                M.best_delay_sec_start(q,s,l) = M.delay_sec_start(b);
                M.best_shape(q,s,l) = M.shape(c);
                M.best_delay_sec_median(q,s,l) = modelwin_convert_delay(M.best_intper_sec(q,s), M.best_delay_sec_start(q,s), M.best_shape(q,s), 'median');
                
                % get prediction
                if (l == 1 || ~I.skipnullpreds) && (s == 1 || ~I.skipsmppreds)
                    for k = valid_seg_durs
                        [winpow, ~, overlap] = win_power_ratio(L.unique_segs(k)/1000, I.distr, ...
                            M.best_intper_sec(q,s), M.best_delay_sec_start(q,s), ...
                            'shape', M.best_shape(q,s), 'tsec', L.lag_t, ...
                            'rampwin', I.rampwin, 'rampdur', I.rampdur, ...
                            'intervalmass', I.intervalmass, ...
                            'intervaltype', I.intervaltype, ...
                            'delaypoint', 'start');
                        if I.winpowratio
                            predictor_notrans = winpow;
                        else
                            predictor_notrans = overlap;
                        end
                        predictor_notrans(predictor_notrans<0) = 0;
                        predictor = tranpred(predictor_notrans);
                        
                        M.diff_context_bestpred(:,k,q,s) = bsxfun(@times, predictor, M.same_context(:,k,q,s));
                    end
                end
            end
        end
    end
    
    % compute significance
    if I.nullsmps>0
        % -> null sample x channel by sample
        X = permute(M.best_loss, [3, 1, 2]);
        
        % average across splits
        if isfield(L, 'splitsmps')
            X = mean(X(:,:,L.splitsmps),3);
        end
        
        % significance
        M.logP_gaussfit = sig_via_null_gaussfit(-X(1,:,:), -X(2:end,:,:));
        M.logP_counts = sig_via_null_gaussfit(-X(1,:,:), -X(2:end,:,:));
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
    
    n_smps = size(M.diff_context,4);
    n_channels = length(L.channels);
    for q = 1:n_channels
        for s = 1:((n_smps-1)*I.plot_extrasmps+1)
            for b = 1:(I.nullsmps*I.plot_nullsmps+1)
                
                chan = L.channels(q);
                
                if strcmp(L.chnames{q}, ['ch' num2str(chan)])
                    chname = L.chnames{q};
                else
                    chname = ['ch' num2str(chan) '-' L.chnames{q}];
                end
                
                if s > 1
                    smpstr = ['_smp' num2str(s-1)];
                else
                    smpstr = '';
                end
                
                if b > 1
                    nullstr = ['_null' num2str(b-1)];
                else
                    nullstr = '';
                end
                
                % plot prediction for best delay, lineplot
                if (s == 1 || ~I.skipsmppreds) && (b == 1 || ~I.skipnullpreds)
                    
                    clf(I.figh);
                    set(I.figh, 'Position', [100, 100, 900, 900]);
                    
                    X = M.diff_context(:,:,q,s,b);
                    corr_range = quantile(X(:), [0.01, 0.99]);
                    clear X;
                    invariance_line = NaN;
                    valid_seg_durs = find(L.n_total_segs(:,1)>0)';
                    for k = valid_seg_durs
                        subplot(4, 2, k);
                        hold on;
                        plot(L.lag_t([1 end]) * 1000, [0,0], 'k--', 'LineWidth', I.linewidth);
                        h1 = plot(L.lag_t * 1000, M.same_context(:,k,q,s,b), 'LineWidth', I.linewidth);
                        h2 = plot(L.lag_t * 1000, M.diff_context(:,k,q,s,b), 'LineWidth', I.linewidth);
                        h3 = plot(L.lag_t * 1000, M.diff_context_bestpred(:,k,q,s,b), 'LineWidth', I.linewidth);
                        plot(L.unique_segs(k)*[1 1], corr_range, 'k--', 'LineWidth', I.linewidth);
                        if ~isnan(invariance_line); plot(I.plot_win * 1000, invariance_line*[1 1], 'k--', 'LineWidth', 2); end
                        xlim(I.plot_win * 1000);
                        ylim(corr_range);
                        xlabel('Time Lag (ms)');
                        ylabel('Pearson Correlation');
                        title(sprintf('Seg: %.0f ms', L.unique_segs(k)))
                        if k == 1
                            legend([h1, h2, h3], 'Same', 'Cross', 'Model');
                        end
                        box off;
                    end
                    fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-prediction-lineplot_plotwin' plot_win_string smpstr nullstr]);
                    % print_wrapper([fname '.png']);
                    % print_wrapper([fname '.pdf']);
                    export_fig([fname '.pdf'], '-pdf', '-transparent');
                    export_fig([fname '.png'], '-png', '-transparent', '-r150');
                    savefig(I.figh, [fname '.fig']);
                    
                    % plot prediction for best delay, image
                    clf(I.figh);
                    set(I.figh, 'Position', [100, 100, 900, 600]);
                    subplot(2,1,1);
                    imagesc(M.diff_context(:,:,q,s,b)', corr_range(2) * [-1, 1]);
                    subplot(2,1,2);
                    imagesc(M.diff_context_bestpred(:,:,q,s,b)', corr_range(2) * [-1, 1]);
                    for i = 1:2
                        subplot(2,1,i);
                        colorbar;
                        colormap(flipud(cbrewer('div', 'RdBu', 128)));
                        set(gca, 'YTick', 1:length(L.unique_segs), 'YTickLabel', L.unique_segs);
                        xticks = get(gca, 'XTick');
                        set(gca, 'XTick', xticks, 'XTickLabel', L.lag_t(xticks)*1000);
                        ylabel('Seg Dur (ms)'); xlabel('Lag (ms)');
                        if i == 2
                            title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q,s,b)*1000, M.best_delay_sec_median(q,s,b)*1000));
                        end
                    end
                    fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-prediction_plotwin' plot_win_string smpstr nullstr]);
                    % print_wrapper([fname '.png']);
                    export_fig([fname '.png'], '-png', '-transparent', '-r150');
                    savefig(I.figh, [fname '.fig']);
                end
                
                % delays to plot
                X_start_delays = M.loss(:,:,M.best_shape(q,s,b)==M.shape,q,s,b);
                if ~strcmp(I.plot_delaystat, 'start')
                    delay_sec_altstat = nan(size(X_start_delays));
                    for l = 1:length(M.intper_sec)
                        delay_sec_altstat(l,:) = modelwin_convert_delay(M.intper_sec(l), M.delay_sec_start, M.best_shape(q,s,b), I.plot_delaystat);
                    end
                    delays_to_plot = I.plot_delay_range(1):I.delay_interval:I.plot_delay_range(2);
                    X_altdelays = nan(length(M.intper_sec), length(delays_to_plot));
                    for l = 1:length(M.intper_sec)
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
                bounds = [minX, quantile(X(:,xi), I.ploterrquant)];
                clear xi zi;
                % bounds = quantile(-X(:), [1-I.ploterrquant, 1]);
                if ~all(isnan(X(:)))
                    clf(I.figh);
                    set(I.figh, 'Position', [100, 100, 600, 600]);
                    imagesc(X, bounds);
                    colormap(cmap);
                    colorbar;
                    yticks = unique(round(linspace(1, length(M.intper_sec), 5)));
                    set(gca, 'YTick', yticks, 'YTickLabel', 1000*M.intper_sec(yticks));
                    % ytick = interp1(log2(M.intper_sec), 1:length(M.intper_sec), log2(M.intper_sec));
                    % set(gca, 'YTick', ytick, 'YTickLabel', M.intper_sec*1000);
                    xticks = unique(round(linspace(1, length(delays_to_plot), 5)));
                    set(gca, 'XTick', xticks, 'XTickLabel', 1000*delays_to_plot(xticks));
                    xlabel('Delay (ms)'); ylabel('Receptive Field Duration (ms)');
                    title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q,s,b)*1000, M.best_delay_sec_median(q,s,b)*1000));
                    set(gca, 'FontSize', 12);
                    fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-error-delay-' I.plot_delaystat smpstr nullstr]);
                    % print_wrapper([fname '.png']);
                    export_fig([fname '.png'], '-png', '-transparent', '-r150');
                    savefig(I.figh, [fname '.fig']);
                    % clear X;
                end
            end
        end
    end
end