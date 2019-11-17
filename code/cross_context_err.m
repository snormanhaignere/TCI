function [L, M, MAT_file] = lag_error_cross_segdur_modelfit_modular(D, t, S, varargin)

% Calculates lag-based correlation to estimate integration periods
%
% D: time x stim x channel x repetition
%
% t: vector of time stamps; length(t) must equal size(D,1)
%
% S: structure with stimulus information
%
% See comments below for optional inputs (all fields of I)
%
% 2019-04-05: Last edited by Sam NH

global root_directory;
root_directory = my_root_directory;
addpath(genpath([root_directory '/export_fig_v3']));
addpath(genpath([root_directory '/general-analysis-code']));

% similarity metric used to compare responses
% options are 'corr' (pearson correlation) 
% or 'nse' (normalized squared error)
I.errfunc = 'square';

% sources to exclude specified as indices into S.sources
I.excludesources = [];

% boundary constraint
% 'none': no constraint
% 'noleft': short seg cannot be on left boundary of long seg
% 'noright': short seg cannot be on right boundary of long seg
% 'nobound': short seg cannot be on either boundary of long seg
% 'yesbound': short seg must be on boundary of the "longer" seg
% 'right': must be on right of boundary
% 'left': must be on left of boundary
% (in practice means we must use same segment duration)
I.boundary = 'none';

% can optionally force longer segment to be greater than a fixed value
I.minlongseg = NaN;

% window to use for computing analysis
% and window to use for plotting
I.lag_win = [0, 1];
I.plot_win = I.lag_win;

% channels for which to perform analysis (default: all)
I.channels = 1:size(D,3);

% whether to plot figures
I.plot_figure = true;

% range of values to plot
% if NaN, chosen based on max value of data
I.plot_range = NaN;

% other plotting parameters
I.winrange = [0, 6]; % in octaves relative to base
I.winbase = 0.03125; % in seconds
I.winhop = 1;
I.winsize = 1;

% modelfit parameters
I.distr = 'gamma';
I.forcecausal = false;
I.bestcausal = true; % select only amongst causal solutions (not relevant for Gausssian)
I.plotcausal = true; % only plot errors for causal solutions
switch I.distr
    case 'gauss'
        I.shape = 1;
    case 'gamma'
        I.shape = logspace(log10(1), log10(10), 10);
    otherwise
        error('No matching distribution');
end

% whether to overwrite or debug
I.keyboard = false;
I.overwrite = false;

% directory to save results to
I.output_directory = pwd;

% can optionally specify alternative directory to save figures to
% default is to use same directory
I.figure_directory = I.output_directory;

% default names for each channel
I.chnames = cell(size(D,3));
for i = 1:size(D,3)
    I.chnames{i} = ['ch' num2str(i)];
end

% if false, no analysis is performed, but MAT file is returned
I.run = true;

% overwrite defaults with optional inputs
[I, C, C_value] = parse_optInputs_keyvalue(varargin, I);

% if figure directory not specified, use output directory
% need to re-assign if output directory changed
if ~optInputs(varargin, 'figure_directory') && C.output_directory
    I.figure_directory = I.output_directory;
end

% set plotting window equal to analysis/lag window
% need to re-assign if lag window has changed
if ~optInputs(varargin, 'plot_win') && C.lag_win
    I.plot_win = I.lag_win;
end

% enter debugging mode
if I.keyboard
    keyboard
end

% parameter and subject string identifying this analysis
param_string = optInputs_to_string(I, C_value, ...
    {'errfunc','boundary','lag_win', 'distr'}, ...
    {'plot_win', 'plot_figure', 'plot_range', ...
    'keyboard', 'overwrite', 'output_directory', 'figure_directory', ...
    'chnames', 'run', 'winrange', 'winbase', 'winhop', 'winsize'});

% % parameter with optional arguments
% param_string = [...
%     struct2string(I, 'include_fields', {'errfunc','boundary','lag_win'}), ...
%     '_' struct2string(C_value, 'include_fields', {'minlongseg','channels','excludesources'})];
% if iscell(param_string) || length(param_string) > 120
%     param_string = optInputs_to_string(I, C_value, ...
%         {'errfunc','boundary','lag_win'}, ...
%         {'plot_win', 'plot_figure', 'plot_range', ...
%         'keyboard', 'overwrite', 'output_directory', 'figure_directory', ...
%         'chnames', 'run', 'winrange', 'winbase', 'winhop', 'winsize'});
% else
%     if param_string(end)=='_'; param_string = param_string(1:end-1); end
% end
        
% directory and MAT file to save results to
MAT_file = mkpdir([I.output_directory '/' param_string '/lagcorr.mat']);

% can stop and just return the MAT_file
% (useful to check if it exists)
if ~I.run || isempty(I.channels)
    % save additional parameters
    clear L;
    L.channels = I.channels;
    L.param_string = param_string;
    L.chnames = I.chnames;
    L.boundary = I.boundary;
    L.figure_directory = [I.figure_directory '/' param_string];
    L.output_directory = [I.output_directory '/' param_string];
    L.sr = 1/(t(2)-t(1));
    L.lag_t = I.lag_win(1):1/L.sr:I.lag_win(2);
    M = [];
    return;
end

% similarity function used to compare responses
switch I.errfunc
    case 'square'
        errfunc = @(a,b)(a-b).^2;
    case 'abs'
        errfunc = @(a,b)abs(a-b);
    otherwise
        error('No matching similarity function')
end

if ~exist(MAT_file, 'file') || I.overwrite
    
    % window over which to compute lags
    L.sr = 1/(t(2)-t(1));
    L.lag_t = I.lag_win(1):1/L.sr:I.lag_win(2);
    n_lags = length(L.lag_t);
    
    % total number of channels to test
    n_channels = length(I.channels);
    
    % remove dimensions with only NaNs
    xi = any(any(any(~isnan(D),1),2),3);
    D_noNaN = D(:,:,I.channels,xi);
    n_reps = sum(xi);
    clear xi;
    
    %% Segment duration ranges for plotting
    
    segdur_centers = I.winbase * 2.^(I.winrange(1):I.winhop:I.winrange(2));
    n_win = length(segdur_centers);
    segdur_ranges = nan(n_win,2);
    for i = 1:n_win
        segdur_ranges(i,:) = segdur_centers(i) * 2.^(I.winsize*[-1 1]/2);
    end
    
    %% Model fit parameters
    
    M.distr = I.distr;
    M.intper_sec = logspace(log10(min(S.flatseginfo.durs)), log10(max(S.flatseginfo.durs)), 20);
    M.delay_smp = 1:0.5*L.sr;
    M.shape = I.shape;
    M.err = zeros(length(M.intper_sec), length(M.delay_smp), length(M.shape), n_channels);
        
    %% Create difference score
    
    max_n_pairs = 100000;
    L.same_context = zeros(n_lags, n_win, n_channels); % response pairs;
    L.diff_context = zeros(n_lags, n_win, n_channels); % response pairs;
    L.ncomparisons = zeros(1, max_n_pairs);
    L.n_segs_per_win = zeros(1, n_win);
    L.segdurs = nan(1, max_n_pairs);
    L.pairs = [];
    pair_count = 0;
    
    % choose which target segs to look for
    n_segs = length(S.flatseginfo.source_onset_times);
    target_segs = 1:n_segs;
    if ~isempty(I.excludesources)
        xi = ismember(S.flatseginfo.source_ids, I.excludesources);
        target_segs = target_segs(~xi);
        clear xi;
    end
    
    % saved predictors
    used_predictors = {};
    used_segdurs = [];
    for i = 1:length(target_segs)
        
        targ_seg = target_segs(i);
        
        if mod(i,10)==0
            fprintf('targ seg %d of %d\n', targ_seg, length(S.flatseginfo.source_onset_times));
            drawnow;
        end
        
        % find other segments where the target segment was present but context was different
        same_scramstim = S.flatseginfo.scrambled_ids == S.flatseginfo.scrambled_ids(targ_seg);
        same_source = S.flatseginfo.source_ids == S.flatseginfo.source_ids(targ_seg);
        contains_seg = ge_tol(S.flatseginfo.source_onset_times(targ_seg), S.flatseginfo.source_onset_times) & le_tol(S.flatseginfo.source_offset_times(targ_seg), S.flatseginfo.source_offset_times);
        all_valid_segments = ~same_scramstim & same_source & contains_seg;
        
        % optionally remove segments that aren't long enough
        if ~isnan(I.minlongseg)
            all_valid_segments = all_valid_segments & ge_tol(S.flatseginfo.durs, I.minlongseg);
        end
        
        % add additional constraints
        switch I.boundary
            case 'none'
            case 'left'
                on_left_bound = eq_tol(S.flatseginfo.source_onset_times(targ_seg), S.flatseginfo.source_onset_times);
                all_valid_segments = all_valid_segments & on_left_bound;
            case 'right'
                on_right_bound = eq_tol(S.flatseginfo.source_offset_times(targ_seg), S.flatseginfo.source_offset_times);
                all_valid_segments = all_valid_segments & on_right_bound;
            case 'noleft'
                on_left_bound = eq_tol(S.flatseginfo.source_onset_times(targ_seg), S.flatseginfo.source_onset_times);
                all_valid_segments = all_valid_segments & ~on_left_bound;
            case 'noright'
                on_right_bound = eq_tol(S.flatseginfo.source_offset_times(targ_seg), S.flatseginfo.source_offset_times);
                all_valid_segments = all_valid_segments & ~on_right_bound;
            case 'nobound'
                on_left_bound = eq_tol(S.flatseginfo.source_onset_times(targ_seg), S.flatseginfo.source_onset_times);
                on_right_bound = eq_tol(S.flatseginfo.source_offset_times(targ_seg), S.flatseginfo.source_offset_times);
                all_valid_segments = all_valid_segments & ~on_right_bound & ~on_left_bound;
            case 'yesbound'
                on_left_bound = eq_tol(S.flatseginfo.source_onset_times(targ_seg), S.flatseginfo.source_onset_times);
                on_right_bound = eq_tol(S.flatseginfo.source_offset_times(targ_seg), S.flatseginfo.source_offset_times);
                all_valid_segments = all_valid_segments & on_right_bound & on_left_bound;
            otherwise
                error('No matching boundary constraint');
        end
        
        % segments matching all criteria
        matching_seg = find(all_valid_segments);
        
        for j = 1:length(matching_seg)
            
            if isempty(L.pairs) || ~any(targ_seg==L.pairs(:,1) & matching_seg(j)==L.pairs(:,2)) % test if already used this pair
                
                % update pair count
                pair_count = pair_count + 1;
                
                % save the duration of the smaller segment
                L.segdurs(pair_count) = S.flatseginfo.durs(targ_seg);
                
                % interpolate timecourse surrounding target segment
                % time x stim x channel x rep
                % -> time x channel x rep
                try
                    targ_scram_id = S.flatseginfo.scrambled_ids(targ_seg);
                    targ_onset = S.flatseginfo.scrambled_onset_times(targ_seg);
                    X = squeeze_dims(D_noNaN(:,targ_scram_id,:,:),2);
                    Xdims = size(X);
                    X = reshape(X, Xdims(1), prod(Xdims(2:end)));
                    Y_targ = reshape(interp1( t, X, L.lag_t' + targ_onset), [n_lags, Xdims(2:end)]);
                    clear X Xdims targ_scram_id targ_onset;
                    
                    % interpolate timecourse surrounding matching segment
                    % interpolate matching seg
                    % time x stim x channel x rep
                    % -> time x channel x rep
                    match_scram_id = S.flatseginfo.scrambled_ids(matching_seg(j));
                    targ_onset_relative_to_match_onset = S.flatseginfo.source_onset_times(targ_seg) - S.flatseginfo.source_onset_times(matching_seg(j));
                    assert(ge_tol(targ_onset_relative_to_match_onset, 0));
                    match_onset = S.flatseginfo.scrambled_onset_times(matching_seg(j)) + targ_onset_relative_to_match_onset;
                    % match_onset = S.flatseginfo.scrambled_onset_times(matching_seg(j));
                    X = squeeze_dims(D_noNaN(:,match_scram_id,:,:),2);
                    Xdims = size(X);
                    X = reshape(X, Xdims(1), prod(Xdims(2:end)));
                    Y_match = reshape(interp1( t, X, L.lag_t' + match_onset), [n_lags, Xdims(2:end)]);
                    clear X Xdims match_scram_id match_onset;
                catch
                    keyboard
                end
                
                same_context = zeros(n_lags, n_channels);
                diff_context = zeros(n_lags, n_channels);
                
                % compute error
                for k = 1:n_reps
                    for l = k+1:n_reps
                        L.ncomparisons(pair_count) = L.ncomparisons(pair_count) + 1;
                        diff_context = diff_context + errfunc(Y_targ(:,:,k), Y_match(:,:,l))/2 + errfunc(Y_targ(:,:,l), Y_match(:,:,k))/2;
                        same_context = same_context + errfunc(Y_targ(:,:,k), Y_targ(:,:,l))/2 + errfunc(Y_match(:,:,k),Y_match(:,:,l))/2;
                    end
                end
                clear Y_targ Y_match;
                                
                % add pairs to library
                L.pairs = [L.pairs; targ_seg, matching_seg(j)];
                L.pairs = [L.pairs; matching_seg(j), targ_seg];
                
                % add this value
                xi = ge_tol(L.segdurs(pair_count),segdur_ranges(:,1)) & le_tol(L.segdurs(pair_count),segdur_ranges(:,2));
                assert(sum(xi)<2);
                if sum(xi>0)
                    L.diff_context(:,xi,:) = squeeze_dims(L.diff_context(:,xi,:),2) + diff_context;
                    L.same_context(:,xi,:) = squeeze_dims(L.same_context(:,xi,:),2) + same_context;
                    L.n_segs_per_win(xi) = L.n_segs_per_win(xi) + 1;
                end
                
                % difference vector to predict with model
                Y = diff_context-same_context;
                
                % load predcitor if already computed
                predictor_found = false;
                if ~isempty(used_segdurs)
                    xi = abs(L.segdurs(pair_count)-used_segdurs)<1e-3;
                    if any(xi)
                        assert(sum(xi)==1);
                        predictor = used_predictors{xi};
                        predictor_found = true;
                    end
                end
                
                % otherwise compute
                if ~predictor_found
                    [predictor, t_sec, M] = model_predictors(...
                        L.segdurs(pair_count), M, [-L.lag_t(end), L.lag_t(end)], L.sr, ...
                        'forcecausal', I.forcecausal);
                    predictor = predictor(t_sec >= 0,:,:,:);
                    predictor = predictor(1:n_lags,:,:,:);
                    predictor = 1-predictor;
                    used_predictors = [used_predictors, {predictor}]; %#ok<AGROW>
                    used_segdurs = [used_segdurs; L.segdurs(pair_count)]; %#ok<AGROW>
                end
                
                for m = 1:length(M.shape)
                    for n = 1:length(M.intper_sec)
                        % fprintf('shape: %.2f, intper: %.2f ms\n', M.shape(m), M.intper_sec(i)*1000);
                        for o = 1:length(M.delay_smp)
                            X = predictor(:,n,o,m);
                            if any(any(isnan(Y)))
                                for q = 1:n_channels
                                    xi = ~isnan(Y(:,q));
                                    Yh = X(xi)*(pinv(X(xi))*Y(xi,q));
                                    M.err(n,o,m,q) = M.err(n,o,m,q) + mean((Yh-Y(xi,q)).^2,1);
                                end
                            else
                                Yh = predictor(:,n,o,m)*(pinv(predictor(:,n,o,m))*Y);
                                M.err(n,o,m,:) = M.err(n,o,m,:) + reshape(mean((Yh-Y).^2,1),[1,1,1,n_channels]);
                            end
                        end
                    end
                end
            end
        end
    end
    
    % remove trailing zeros/NaNs
    L.ncomparisons = L.ncomparisons(1:pair_count);
    L.segdurs = L.segdurs(1:pair_count);
    
    % find best prediction
    M.best_intper_sec = nan(n_channels,1);
    M.best_delay_smp = nan(n_channels,1);
    M.best_shape_smp = nan(n_channels,1);
    for q = 1:n_channels
        X = M.err(:,:,:,q);
        if any(M.causal_win(:)) && I.bestcausal
            X(~M.causal_win) = inf;
        end
        [~,xi] = min(X(:));
        [a,b,c] = ind2sub(size(X), xi);
        M.best_intper_sec(q) = M.intper_sec(a);
        M.best_delay_smp(q) = M.delay_smp(b);
        M.best_shape(q) = M.shape(c);
        clear X;
    end
    

    % save
    save(MAT_file, 'L', 'M', '-v7.3');
    
else
    
    load(MAT_file, 'L', 'M');
    
end


% save additional parameters
L.channels = I.channels;
M.channels = L.channels;
L.param_string = param_string;
L.chnames = I.chnames;
L.boundary = I.boundary;
L.figure_directory = [I.figure_directory '/' param_string];
L.output_directory = [I.output_directory '/' param_string];

%% Plotting
% center of the windows
if I.plot_figure
    
    figh = figure;
            
    segdur_centers = I.winbase * 2.^(I.winrange(1):I.winhop:I.winrange(2));
    n_win = length(segdur_centers);
    n_channels = size(L.same_context,3);
    
    for q = 1:n_channels
        
        chan = L.channels(q);
        if strcmp(L.chnames{chan}, ['ch' num2str(chan)])
            chname = L.chnames{chan};
        else
            chname = ['ch' num2str(chan) '-' L.chnames{chan}];
        end
        
        % difference
        n_lags = length(L.lag_t);
        Z = nan(n_lags, n_win);
        Zp = nan(n_lags, n_win);
        for i = 1:n_win
            Z(:,i) = (L.diff_context(:,i,q) - L.same_context(:,i,q))/L.n_segs_per_win(i);
            
            [overlap,~,t_sec] = winoverlap(segdur_centers(i), ...
                I.distr, M.best_intper_sec(q), M.best_delay_smp(q)/L.sr, 'shape', M.best_shape(q), ...
                'win', [-L.lag_t(end), L.lag_t(end)], 'sr', L.sr, 'forcecausal', I.forcecausal, 'plot', false);
            predictor = overlap(t_sec >= 0);
            predictor = predictor(1:n_lags);
            predictor = predictor(:);
            predictor = 1-predictor;
            xi = ~isnan(Z(:,i));
            Zp(xi,i) = predictor(xi)*pinv(predictor(xi))*Z(xi,i);
        end
                
        %%
        
        clf(figh);
        set(figh, 'Position', [100, 100, 900, 600]);
        n_rows = ceil(sqrt(n_win));
        n_cols = ceil(n_win/n_rows);
        ybounds = quantile(Z(:), [0, 0.99]);
        ybounds(1) = min(ybounds(1), 0);
%         ybounds = [-0.001, 0.01]
        for i = 1:n_win
            subplot(n_rows, n_cols, i);
            plot(L.lag_t*1000, Z(:,i), 'LineWidth', 2);
            hold on;
            plot(L.lag_t*1000, Zp(:,i), 'r-', 'LineWidth', 2);
            plot(L.lag_t*1000, zeros(size(L.lag_t)), 'k-', 'LineWidth', 2);
            plot(segdur_centers(i)*[1 1]*1000, ybounds, 'r--', 'LineWidth', 2);
            ylim(ybounds);
            xlim(I.plot_win*1000);
            title(sprintf('%.0f ms', segdur_centers(i)*1000));
            xlabel('Time (ms)');
            ylabel('Error');
        end
        
        export_fig(mkpdir([L.figure_directory '/' chname ...
            '-win' num2str(I.plot_win(1)) '-' num2str(I.plot_win(2))]), '-pdf', '-transparent');
        
        
        %%
        
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
            clf(figh);
            set(figh, 'Position', [100, 100, 600, 600]);
            
            imagesc(-X, bounds);
            colormap(parula(128));
            colorbar;
            ytick = interp1(log2(M.intper_sec), 1:length(M.intper_sec), log2(M.intper_sec));
            set(gca, 'YTick', ytick, 'YTickLabel', M.intper_sec*1000);
            xticks = get(gca, 'XTick');
            set(gca, 'XTickLabel', 1000*M.delay_smp(xticks)/L.sr);
            xlabel('Delay (ms)'); ylabel('Receptive Field Duration (ms)');
            title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q)*1000, M.best_delay_smp(q)/L.sr*1000));
            set(gca, 'FontSize', 12);
            fname = [L.figure_directory '/' chname '-model-error_' causal_string];
            print_wrapper([fname '.png']);
            % export_fig([fname '.png'], '-png', '-transparent', '-r100');
            clear X;
        end
    end
end
