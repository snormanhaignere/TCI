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
I.corrfac = 0;

% whether to transform the weighting applied to segments
I.tranweightnsegs = 'none'; % applied to total segs

% whether to transform the predictions
I.tranpred = 'none';

% whether to create null samples with phase scrambling
I.nullsmps = 0;

% skip calculation of the best predictors for null samples
I.skipnullpreds = true;

% skip calculation of best predictors for extra samples
I.skipsmppreds = ~isfield(L, 'splits_diff_context');

% whether to skip cross-validated predictions
I.skipcvpreds = false;

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

% smoothing window for plotting
I.plot_smoothwin = 0;

% line width for plotting
I.linewidth = 2;

% quantile of error map to plot
I.ploterrquant = 0.3;

% whether to plot splits or resamples
I.plot_splits = isfield(L, 'splits_diff_context');
I.plot_bstrap = false;
I.plot_jack = false;
I.plot_nullsmps = false;

% directory to save output to
I.output_directory = L.output_directory;

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
    'run', 'figh', 'keyboard', 'plot_figure', 'plot_win', 'plot_smoothwin', 'linewidth', ...
    'overwrite', 'ploterrquant', 'plot_splits', 'plot_delaystat', 'plot_delay_range', ...
    'output_directory', 'plot_bstrap', 'plot_jack'};
param_string_modelfit = optInputs_to_string(I, C_value, always_include, always_exclude);

% file to save results to
% can optionally and just return this file
MAT_file = mkpdir([I.output_directory '/model_fit_' param_string_modelfit '.mat']);
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

[n_lags, n_seg, n_channels, ~] = size(L.same_context);
assert(n_channels == length(L.channels));
assert(length(L.unique_segs)==n_seg);
assert(n_lags == length(L.lag_t))

%% Model fit

if ~exist(MAT_file, 'file') || I.overwrite
    
    clear M;
    ResetRandStream2(I.randseed);
    
    % assign relevant structures from L
    M.diff_context = L.diff_context;
    M.same_context = L.same_context;
    M.same_context_err = L.same_context_err;
    M.n_total_segs = L.n_total_segs;
    
    % add splits to this same structure
    % since were going to do the same thing
    % for these splits, and its more efficient
    % to do everything once
    % data from splits are separated out at the end
    if isfield(L, 'splits_diff_context')
        [~, ~, ~, n_partitions, n_splits] = size(L.splits_diff_context);
        split_smps = (1:(n_partitions*n_splits))+1;
        d = [n_lags, n_seg, n_channels, n_partitions*n_splits];
        M.diff_context = cat(4, M.diff_context, reshape(L.splits_diff_context, d));
        M.same_context = cat(4, M.same_context, reshape(L.splits_same_context, d));
        M.same_context_err = cat(4, M.same_context_err, reshape(L.splits_same_context_err, d));
        M.n_total_segs = cat(2, M.n_total_segs, reshape(L.splits_n_total_segs, [n_seg, n_partitions*n_splits]));
        clear d;
    end
    
    % add bootstraps to the same structure
    if isfield(L, 'bstrap_diff_context')
        n_bstrap = size(L.bstrap_diff_context,4);
        bstrap_smps = (1:n_bstrap)+size(M.diff_context,4);
        M.diff_context = cat(4, M.diff_context, L.bstrap_diff_context);
        M.same_context = cat(4, M.same_context, L.bstrap_same_context);
        M.same_context_err = cat(4, M.same_context_err, L.bstrap_same_context_err);
        M.n_total_segs = cat(2, M.n_total_segs, L.bstrap_n_total_segs);
    end
    
    % add jack-knife to the same structure
    if isfield(L, 'jack_diff_context')
        n_jack = size(L.jack_diff_context,4);
        jack_smps = (1:n_jack)+size(M.diff_context,4);
        M.diff_context = cat(4, M.diff_context, L.jack_diff_context);
        M.same_context = cat(4, M.same_context, L.jack_same_context);
        M.same_context_err = cat(4, M.same_context_err, L.jack_same_context_err);
        M.n_total_segs = cat(2, M.n_total_segs, L.jack_n_total_segs);
    end
    n_smps = size(M.diff_context,4);
    
    % optionally create null samples of cross-context correlation via phase
    % scrambling, leave same context correlations unaltered by simply
    % repeating them over the null sample dimension
    % again we append everything, now to a new fifth dimension
    % to save computation time
    if I.nullsmps>0
        diff_context_null = nan([n_lags, n_seg, n_channels, n_smps, I.nullsmps]);
        for i = 1:I.nullsmps
            diff_context_null(:,:,:,:,i) = real(phasescram(M.diff_context));
        end
        
        %         keyboard
        %         l = 4;
        %         d = M.diff_context(:,l,2,13,1);
        %         dn = squeeze(diff_context_null(:,l,2,13,2:3));
        %         figure;
        %         plot(L.lag_t, dn)
        %         hold on;
        %         plot(L.lag_t, d, 'k-', 'LineWidth', 3)
        %         sum(d.^2)
        %         sum(dn.^2)
        
        assert(size(M.diff_context,5)==1);
        M.diff_context = cat(5, M.diff_context, diff_context_null);
        M.same_context = repmat(M.same_context, [1, 1, 1, 1, I.nullsmps+1]);
        M.same_context_err = repmat(M.same_context_err, [1, 1, 1, 1, I.nullsmps+1]);
        clear diff_context_null;
    end
    
    % valid segment durations
    valid_seg_durs = find(M.n_total_segs(:,1)>0)';
    n_valid_segdurs = length(valid_seg_durs);
    
    % segment dependent weights
    % segdur x data sample
    W_segs = tranweightnsegsfn(M.n_total_segs);
    
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
                    'delaypoint', 'start',  'corrfac', I.corrfac);
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
                    %                     keyboard;
                    %
                    %                     %%
                    %
                    %                     figure;
                    %                     imagesc(squeeze(diff_context_pred(:,2,13,:)*1e6));
                    %                     colorbar;
                    %
                    %                     %%
                    %
                    %                     x = squeeze(sum(diff_context_format(1:100,2,13,:).^2))
                    %
                    %
                    %                     %%
                    %
                    %                     e = x-y;
                    %
                    %
                    %                     %%
                    %
                    %                     plot(W_total_format(:,2,13,1))
                    %
                    %                     %%
                    %
                    %                     figure;
                    %                     imagesc(squeeze(E(1,2,:,:)))
                    %
                    %                     %%
                    M.loss(i,j,m,:,:,:) = lossfn(diff_context_format, diff_context_pred, W_total_format);
                end
            end
        end
    end
    
    % find best prediction
    M.diff_context_bestpred = nan(n_lags, n_seg, n_channels, (n_smps-1)*(~I.skipsmppreds)+1, I.nullsmps*(~I.skipnullpreds)+1);
    M.best_intper_sec = nan(n_channels, n_smps, I.nullsmps+1);
    M.best_delay_sec_start = nan(n_channels, n_smps, I.nullsmps+1);
    M.best_delay_sec_median = nan(n_channels, n_smps, I.nullsmps+1);
    M.best_shape = nan(n_channels, n_smps, I.nullsmps+1);
    M.best_loss = nan(n_channels, n_smps, I.nullsmps+1);
    for q = 1:n_channels
        for s = 1:n_smps
            for l = 1:I.nullsmps+1
                X = M.loss(:,:,:,q,s,l);
                [~,xi] = min(X(:));
                [a,b,c] = ind2sub(size(X), xi);
                clear X;
                M.best_loss(q,s,l) = M.loss(a,b,c,q,s,l);
                M.best_intper_sec(q,s,l) = M.intper_sec(a);
                M.best_delay_sec_start(q,s,l) = M.delay_sec_start(b);
                M.best_shape(q,s,l) = M.shape(c);
                M.best_delay_sec_median(q,s,l) = modelwin_convert_delay(M.best_intper_sec(q,s,l), ...
                    M.best_delay_sec_start(q,s,l), M.best_shape(q,s,l), 'median', ...
                    'intervaltype', I.intervaltype, 'intervalmass', I.intervalmass);
                clear a b c;
                
                % get prediction
                if (l == 1 || ~I.skipnullpreds) && (s == 1 || ~I.skipsmppreds)
                    for k = valid_seg_durs
                        [winpow, ~, overlap] = win_power_ratio(L.unique_segs(k)/1000, I.distr, ...
                            M.best_intper_sec(q,s,l), M.best_delay_sec_start(q,s,l), ...
                            'shape', M.best_shape(q,s,l), 'tsec', L.lag_t, ...
                            'rampwin', I.rampwin, 'rampdur', I.rampdur, ...
                            'intervalmass', I.intervalmass, ...
                            'intervaltype', I.intervaltype, ...
                            'delaypoint', 'start', 'corrfac', I.corrfac);
                        if I.winpowratio
                            predictor_notrans = winpow;
                        else
                            predictor_notrans = overlap;
                        end
                        predictor_notrans(predictor_notrans<0) = 0;
                        predictor = tranpred(predictor_notrans);
                        
                        M.diff_context_bestpred(:,k,q,s,l) = bsxfun(@times, predictor, M.same_context(:,k,q,s,l));
                    end
                end
            end
        end
    end
    
    % separate out null samples into a different field for clarity
    if I.nullsmps>0
        fa = {...
            'loss', 6; ...
            'diff_context', 5; 'same_context', 5; 'same_context_err', 5; ...
            'best_intper_sec', 3; 'best_delay_sec_start', 3; 'best_delay_sec_median', 3; 'best_shape', 3; 'best_loss', 3; ...
            };
        if ~I.skipnullpreds
            fa = cat(1, fa, {'diff_context_bestpred'; 5});
        end
        for j = 1:length(fa)
            f = fa{j,1};
            dim = fa{j,2};
            M = split_dim_into_fields(M, f, dim, {1, (1:I.nullsmps)+1}, [{f}, strcat('null_', {f})]);
        end
    end
    
    % compute significance using the null samples
    % -> null sample x channel by sample
    if I.nullsmps>0
        M.logP_gaussfit = sig_via_null_gaussfit(-M.best_loss, -permute(M.null_best_loss,[3,1,2]));
        M.logP_counts = sig_via_null_counts(-M.best_loss, -permute(M.null_best_loss,[3,1,2]));
    end
    
    % separate out partitions/splits, bootstraps, and jack-knifes
    if n_smps > 1
        smp_groups = {};
        prefixes = {};
        shapes = {};
        if isfield(L, 'splits_diff_context')
            smp_groups = [smp_groups, {split_smps}];
            prefixes = [prefixes, {'splits_'}];
            shapes = [shapes, {[n_partitions, n_splits]}];
        end
        if isfield(L, 'bstrap_diff_context')
            smp_groups = [smp_groups, {bstrap_smps}];
            prefixes = [prefixes, {'bstrap_'}];
            shapes = [shapes, {n_bstrap}];
        end
        if isfield(L, 'jack_diff_context')
            smp_groups = [smp_groups, {jack_smps}];
            prefixes = [prefixes, {'jack_'}];
            shapes = [shapes, {n_jack}];
        end
        
        % fields and dimensions to separate out
        fa = {...
            'loss', 5; ...
            'diff_context', 4; 'same_context', 4; 'same_context_err', 4; ...
            'best_intper_sec', 2; 'best_delay_sec_start', 2; 'best_delay_sec_median', 2; 'best_shape', 2; 'best_loss', 2; ...
            };

        if I.nullsmps>0
            n_f = length(fa);
            for i = 1:n_f
                fa = cat(1, fa, {['null_' fa{i,1}], fa{i, 2}});
            end
            fa = cat(1, fa, {'logP_gaussfit', 2; 'logP_counts', 2});
        end
        if ~I.skipsmppreds
            fa = cat(1, fa, {'diff_context_bestpred', 4});
            if ~I.skipnullpreds
                fa = cat(1, fa, {'null_diff_context_bestpred', 4});
            end
        end
        for j = 1:length(fa)
            f = fa{j,1};
            dim = fa{j,2};
            M = split_dim_into_fields(M, f, dim, [{1}, smp_groups], [{f}, strcat(prefixes, f)], ...
                'shapes', [1, shapes], 'squeeze', [true, false(1, length(shapes))]);
        end

    end
    
    % cross-validated loss and significance using independent partitions
    if isfield(L, 'splits_diff_context')
        if ~I.skipcvpreds
            M.cv_diff_context_bestpred = nan(n_lags, n_seg, n_channels, n_partitions, n_splits);
        end
        M.cv_loss = nan(n_channels, n_partitions, n_splits);
        if I.nullsmps>0
            M.null_cv_loss = nan(n_channels, n_partitions, n_splits, I.nullsmps);
        end
        for q = 1:n_channels
            for s = 1:n_splits
                for p = 1:n_partitions
                    X = mean(M.splits_loss(:,:,:,q,setdiff(1:n_partitions,p),s),5);
                    [~,xi] = min(X(:));
                    [a,b,c] = ind2sub(size(X), xi);
                    M.cv_loss(q, p, s) = M.splits_loss(a,b,c,q,p,s);
                    if I.nullsmps>0
                        M.null_cv_loss(q, p, s, :) = M.splits_null_loss(a,b,c,q,p,s,:);
                    end
                    
                    if ~I.skipcvpreds
                        for k = valid_seg_durs
                            [winpow, ~, overlap] = win_power_ratio(L.unique_segs(k)/1000, I.distr, ...
                                M.intper_sec(a), M.delay_sec_start(b), ...
                                'shape', M.shape(c), 'tsec', L.lag_t, ...
                                'rampwin', I.rampwin, 'rampdur', I.rampdur, ...
                                'intervalmass', I.intervalmass, ...
                                'intervaltype', I.intervaltype, ...
                                'delaypoint', 'start', 'corrfac', I.corrfac);
                            if I.winpowratio
                                predictor_notrans = winpow;
                            else
                                predictor_notrans = overlap;
                            end
                            predictor_notrans(predictor_notrans<0) = 0;
                            predictor = tranpred(predictor_notrans);
                            
                            M.cv_diff_context_bestpred(:,k,q,p,s) = bsxfun(@times, predictor, M.splits_same_context(:,k,q,p,s));
                        end
                    end
                    clear a b c X;
                    
                end
            end
        end
        
        % average across partitions and splits
        M.av_cv_loss = mean(mean(M.cv_loss, 2),3);
        if I.nullsmps>0
            M.av_null_cv_loss = squeeze_dims(mean(mean(M.null_cv_loss, 2), 3), [2,3]);
        end
        
        % compute significance
        if I.nullsmps>0
            M.cv_logP_gaussfit = sig_via_null_gaussfit(-M.av_cv_loss, -M.av_null_cv_loss');
            M.cv_logP_counts = sig_via_null_counts(-M.av_cv_loss, -M.av_null_cv_loss');
        end
    end
    
    
    M.channels = L.channels;
    
    save(MAT_file, 'M', '-v7.3');
    
else
    
    load(MAT_file, 'M')
    
end

M.param_string_modelfit = param_string_modelfit;
M.intervaltype = I.intervaltype;
M.intervalmass = I.intervalmass;

%% Plot the results

if I.plot_figure
    
    % create a figure handle if not supplied
    if isempty(I.figh)
        I.figh = figure;
    end
    
    n_channels = length(L.channels);
    for q = 1:n_channels
        
        % channel name
        chan = L.channels(q);
        if strcmp(L.chnames{q}, ['ch' num2str(chan)])
            chname = L.chnames{q};
        else
            chname = ['ch' num2str(chan) '-' L.chnames{q}];
        end
        
        % plot
        aux_args = {M.intper_sec, M.delay_sec_start, L.unique_segs, L.lag_t, I.intervaltype, I.intervalmass, ...
            I.plot_win, I.plot_smoothwin, I.plot_delaystat, I.plot_delay_range, I.ploterrquant, I.linewidth, I.figh};
        fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname]);
        plot_modelfit(M.diff_context(:,:,q), M.same_context(:,:,q), M.diff_context_bestpred(:,:,q), ...
            M.loss(:,:,M.best_shape(q)==M.shape, q), M.best_intper_sec(q), M.best_delay_sec_median(q), ...
            M.best_shape(q), fname, aux_args{:})
        
        % plot null samples
        if I.nullsmps>0 && I.plot_nullsmps
            for l = 1:I.nullsmps
                fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-nullsmp' num2str(l)]);
                if ~I.skipnullpreds
                    null_diff_context_bestpred = M.null_diff_context_bestpred(:,:,q,l);
                else
                    null_diff_context_bestpred = [];
                end
                plot_modelfit(M.null_diff_context(:,:,q,l), M.null_same_context(:,:,q,l), null_diff_context_bestpred, ...
                    M.null_loss(:,:,M.null_best_shape(q)==M.shape,q,l), M.null_best_intper_sec(q,l), ...
                    M.null_best_delay_sec_median(q,l), M.null_best_shape(q,l), fname, aux_args{:})
            end
        end
        
        % plot splits
        if I.plot_splits
            [~, n_partitions, n_splits] = size(M.splits_best_intper_sec);
            for s = 1:n_splits
                for p = 1:n_partitions
                    if ~I.skipcvpreds
                        cv_diff_context_bestpred = M.cv_diff_context_bestpred(:,:,q,p,s);
                    else
                        cv_diff_context_bestpred = [];
                    end
                    fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' ...
                        chname '-split' num2str(s) '-p' num2str(p)]);
                    plot_modelfit(M.splits_diff_context(:,:,q,p,s), M.splits_same_context(:,:,q,p,s), cv_diff_context_bestpred, ...
                        M.splits_loss(:,:,M.best_shape(q)==M.shape,q,p,s), M.splits_best_intper_sec(q,p,s), ...
                        M.splits_best_delay_sec_median(q,p,s), M.splits_best_shape(q,p,s), fname, aux_args{:})
                    %                     if I.nullsmps>0 && I.plot_nullsmps
                    %                         for l = 1:I.nullsmps
                    %                             if ~I.skipsmppreds && ~I.skipnullpreds
                    %                                 splits_null_diff_context_bestpred = M.splits_null_diff_context_bestpred(:,:,q,p,s,l);
                    %                             else
                    %                                 splits_null_diff_context_bestpred = [];
                    %                             end
                    %                             fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' ...
                    %                                 chname '-split' num2str(s) '-p' num2str(p) '-nullsmp' num2str(l)]);
                    %                             plot_modelfit(M.splits_null_diff_context(:,:,q,p,s,l), M.splits_null_same_context(:,:,q,p,s,l), splits_null_diff_context_bestpred, ...
                    %                                 M.splits_null_loss(:,:,M.best_shape(q)==M.shape,q,p,s,l), M.splits_null_best_intper_sec(q,p,s,l), ...
                    %                                 M.splits_null_best_delay_sec_median(q,p,s,l), M.splits_null_best_shape(q,p,s,l), fname, aux_args{:})
                    %
                    %                         end
                    %                     end
                end
            end
        end
        
        
        
        % plot bootstraps
        if I.plot_bstrap
            n_bstrap = size(M.bstrap_diff_context,4);
            for b = 1:n_bstrap
                if ~I.skipsmppreds
                    bstrap_diff_context_bestpred = M.bstrap_diff_context_bestpred(:,:,q,b);
                else
                    bstrap_diff_context_bestpred = [];
                end
                fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' ...
                    chname '-bstrap' num2str(b)]);
                plot_modelfit(M.bstrap_diff_context(:,:,q,b), M.bstrap_same_context(:,:,q,b), bstrap_diff_context_bestpred, ...
                    M.bstrap_loss(:,:,M.bstrap_best_shape(q,b)==M.shape,q,b), M.bstrap_best_intper_sec(q,b), ...
                    M.bstrap_best_delay_sec_median(q,b), M.bstrap_best_shape(q,b), fname, aux_args{:});
            end
        end

        % plot jack-knife samples
        if I.plot_jack
            n_jack = size(M.jack_diff_context,4);
            for b = 1:n_jack
                if ~I.skipsmppreds
                    jack_diff_context_bestpred = M.jack_diff_context_bestpred(:,:,q,b);
                else
                    jack_diff_context_bestpred = [];
                end
                fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' ...
                    chname '-jack' num2str(b)]);
                plot_modelfit(M.jack_diff_context(:,:,q,b), M.jack_same_context(:,:,q,b), jack_diff_context_bestpred, ...
                    M.jack_loss(:,:,M.jack_best_shape(q,b)==M.shape,q,b), M.jack_best_intper_sec(q,b), ...
                    M.jack_best_delay_sec_median(q,b), M.jack_best_shape(q,b), fname, aux_args{:});
            end
        end
        
        %%
        %
        %         for s = 1:((n_smps-1)*I.plot_extrasmps+1)
        %             for b = 1:(I.nullsmps*I.plot_nullsmps+1)
        %
        %                 chan = L.channels(q);
        %
        %                 if strcmp(L.chnames{q}, ['ch' num2str(chan)])
        %                     chname = L.chnames{q};
        %                 else
        %                     chname = ['ch' num2str(chan) '-' L.chnames{q}];
        %                 end
        %
        %                 if s > 1
        %                     if isfield(L, 'spitsmps') && s == n_smps
        %                         smpstr = '_splitav';
        %                     else
        %                         smpstr = ['_smp' num2str(s-1)];
        %                     end
        %                 else
        %                     smpstr = '';
        %                 end
        %
        %                 if b > 1
        %                     nullstr = ['_null' num2str(b-1)];
        %                 else
        %                     nullstr = '';
        %                 end
        %
        %                 % plot prediction for best delay, lineplot
        %                 if (s == 1 || ~I.skipsmppreds) && (b == 1 || ~I.skipnullpreds)
        %
        %                     clf(I.figh);
        %                     set(I.figh, 'Position', [100, 100, 900, 900]);
        %
        %                     X = M.diff_context(:,:,q,s,b);
        %                     corr_range = quantile(X(:), [0.01, 0.99]);
        %                     clear X;
        %                     invariance_line = NaN;
        %                     valid_seg_durs = find(L.n_total_segs(:,1)>0)';
        %                     for k = valid_seg_durs
        %                         subplot(4, 2, k);
        %                         hold on;
        %                         plot(L.lag_t([1 end]) * 1000, [0,0], 'k--', 'LineWidth', I.linewidth);
        %                         h1 = plot(L.lag_t * 1000, M.same_context(:,k,q,s,b), 'LineWidth', I.linewidth);
        %                         h2 = plot(L.lag_t * 1000, M.diff_context(:,k,q,s,b), 'LineWidth', I.linewidth);
        %                         h3 = plot(L.lag_t * 1000, M.diff_context_bestpred(:,k,q,s,b), 'LineWidth', I.linewidth);
        %                         plot(L.unique_segs(k)*[1 1], corr_range, 'k--', 'LineWidth', I.linewidth);
        %                         if ~isnan(invariance_line); plot(I.plot_win * 1000, invariance_line*[1 1], 'k--', 'LineWidth', 2); end
        %                         xlim(I.plot_win * 1000);
        %                         ylim(corr_range);
        %                         xlabel('Time Lag (ms)');
        %                         ylabel('Pearson Correlation');
        %                         title(sprintf('Seg: %.0f ms', L.unique_segs(k)))
        %                         if k == 1
        %                             legend([h1, h2, h3], 'Same', 'Cross', 'Model');
        %                         end
        %                         box off;
        %                     end
        %                     fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-prediction-lineplot_plotwin' plot_win_string smpstr nullstr]);
        %                     % print_wrapper([fname '.png']);
        %                     % print_wrapper([fname '.pdf']);
        %                     export_fig([fname '.pdf'], '-pdf', '-transparent');
        %                     export_fig([fname '.png'], '-png', '-transparent', '-r150');
        %                     savefig(I.figh, [fname '.fig']);
        %
        %                     % plot prediction for best delay, image
        %                     clf(I.figh);
        %                     set(I.figh, 'Position', [100, 100, 900, 600]);
        %                     subplot(2,1,1);
        %                     imagesc(M.diff_context(:,:,q,s,b)', corr_range(2) * [-1, 1]);
        %                     subplot(2,1,2);
        %                     imagesc(M.diff_context_bestpred(:,:,q,s,b)', corr_range(2) * [-1, 1]);
        %                     for i = 1:2
        %                         subplot(2,1,i);
        %                         colorbar;
        %                         colormap(flipud(cbrewer('div', 'RdBu', 128)));
        %                         set(gca, 'YTick', 1:length(L.unique_segs), 'YTickLabel', L.unique_segs);
        %                         xticks = get(gca, 'XTick');
        %                         set(gca, 'XTick', xticks, 'XTickLabel', L.lag_t(xticks)*1000);
        %                         ylabel('Seg Dur (ms)'); xlabel('Lag (ms)');
        %                         if i == 2
        %                             title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q,s,b)*1000, M.best_delay_sec_median(q,s,b)*1000));
        %                         end
        %                     end
        %                     fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-prediction_plotwin' plot_win_string smpstr nullstr]);
        %                     % print_wrapper([fname '.png']);
        %                     export_fig([fname '.png'], '-png', '-transparent', '-r150');
        %                     savefig(I.figh, [fname '.fig']);
        %                 end
        %
        %                 % delays to plot
        %                 X_start_delays = M.loss(:,:,M.best_shape(q,s,b)==M.shape,q,s,b);
        %                 if ~strcmp(I.plot_delaystat, 'start')
        %                     delay_sec_altstat = nan(size(X_start_delays));
        %                     for l = 1:length(M.intper_sec)
        %                         delay_sec_altstat(l,:) = modelwin_convert_delay(M.intper_sec(l), M.delay_sec_start, M.best_shape(q,s,b), I.plot_delaystat);
        %                     end
        %                     delays_to_plot = I.plot_delay_range(1):I.delay_interval:I.plot_delay_range(2);
        %                     X_altdelays = nan(length(M.intper_sec), length(delays_to_plot));
        %                     for l = 1:length(M.intper_sec)
        %                         xi = delays_to_plot > delay_sec_altstat(l,1) & delays_to_plot < delay_sec_altstat(l,end);
        %                         X_altdelays(l,xi) = interp1(delay_sec_altstat(l,:), X_start_delays(l,:), delays_to_plot(xi));
        %                     end
        %                     X_altdelays(isnan(X_altdelays)) = inf;
        %                     X = X_altdelays;
        %                 else
        %                     X = X_start_delays;
        %                 end
        %
        %                 % plot the error vs. parameters
        %                 cmap = flipud(cbrewer('seq', 'Reds', 128));
        %                 [minX,zi] = min(X(:));
        %                 [~, xi] = ind2sub(size(X), zi);
        %                 bounds = [minX, quantile(X(:,xi), I.ploterrquant)];
        %                 clear xi zi;
        %                 % bounds = quantile(-X(:), [1-I.ploterrquant, 1]);
        %                 if ~all(isnan(X(:)))
        %                     clf(I.figh);
        %                     set(I.figh, 'Position', [100, 100, 600, 600]);
        %                     imagesc(X, bounds);
        %                     colormap(cmap);
        %                     colorbar;
        %                     yticks = unique(round(linspace(1, length(M.intper_sec), 5)));
        %                     set(gca, 'YTick', yticks, 'YTickLabel', 1000*M.intper_sec(yticks));
        %                     % ytick = interp1(log2(M.intper_sec), 1:length(M.intper_sec), log2(M.intper_sec));
        %                     % set(gca, 'YTick', ytick, 'YTickLabel', M.intper_sec*1000);
        %                     xticks = unique(round(linspace(1, length(delays_to_plot), 5)));
        %                     set(gca, 'XTick', xticks, 'XTickLabel', 1000*delays_to_plot(xticks));
        %                     xlabel('Delay (ms)'); ylabel('Receptive Field Duration (ms)');
        %                     title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q,s,b)*1000, M.best_delay_sec_median(q,s,b)*1000));
        %                     set(gca, 'FontSize', 12);
        %                     fname = mkpdir([L.figure_directory '/model-fit-' param_string_modelfit '/' chname '-model-error-delay-' I.plot_delaystat smpstr nullstr]);
        %                     % print_wrapper([fname '.png']);
        %                     export_fig([fname '.png'], '-png', '-transparent', '-r150');
        %                     savefig(I.figh, [fname '.fig']);
        %                     % clear X;
        %                 end
        %             end
        %         end
    end
end