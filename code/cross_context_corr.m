function [L, MAT_file] = cross_context_corr(D, t, S, varargin)

% Calculates cross and same context correlation for different lags
% relatice to segment onset.
%
% D: time x stim x channel x repetition
%
% t: vector of time stamps; length(t) must equal size(D,1)
%
% S: structure with stimulsus information
%
% See comments below for optional inputs (all fields of I)

%% Default parameters

clear I;

% similarity metric used to compare responses
% across contexts
% options are:
% 'corr' (pearson correlation) the default
% 'nse' (normalized squared error)
I.simfunc = 'corr';

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

% window of lags computed
I.lag_win = [0, 2*max(S.segs(:))/1000];

% window used for plotting
I.plot_win = I.lag_win;

% can optionally force the longer of the two segment durations being
% compared to have a fixed value
I.longsegdurs = NaN;

% channels for which to perform analysis (default: all)
I.channels = 1:size(D,3);

% whether to plot figures
I.plot_figure = true;
I.plot_adaptation = true;

% range of values to plot
% if NaN, chosen based on max value of data
I.plot_range = NaN;

% whether to overwrite or debug
I.keyboard = false;
I.overwrite = false;

% can exclude segments that come from particular sources
% the sources are specified in S.sources
% specify the indices of the sources to exclude
I.excludesources = [];

% function to transfrom correlation values
% options are 'none', 'z' (ztransform)
I.trancorr = 'none';

% function used to transform weights
% options are 'none', 'sqrt', 'sqrtm3' (square root minus three)
I.tranweight = 'none';

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

% number of permutations
I.nperms = 0;

% whether to interleave short and long segments
I.interleave = false;

% used for debugging purposes
I.simoffset = false;

% figure handle
I.figh = matlab.ui.Figure.empty;

%% Parse user-specified parameters, create parameter string for this analysis

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

% parameter with optional arguments
always_include = {'boundary','lag_win'};
always_exclude = {...
    'plot_win', 'plot_figure', 'plot_adaptation', 'plot_range', ...
    'keyboard', 'overwrite', 'output_directory', 'figure_directory', ...
    'chnames','run','figh'};
param_string = optInputs_to_string(I, C_value, always_include, always_exclude);

% directory and MAT file to save results to
MAT_file = mkpdir([I.output_directory '/' param_string '/corr_data.mat']);

% can stop and just return the MAT_file
% (useful to check if it exists)
if ~I.run
    clear L;
    L.channels = I.channels;
    L.param_string = param_string;
    L.chnames = I.chnames;
    L.boundary = I.boundary;
    L.figure_directory = [I.figure_directory '/' param_string];
    L.output_directory = [I.output_directory '/' param_string];
    L.rampwin = S.rampwin;
    L.rampdur = S.rampdur;
    L.sr = 1/(t(2)-t(1));
    L.lag_t = I.lag_win(1):1/L.sr:I.lag_win(2);    
    return;
end

% function used to transform correlation
switch I.trancorr
    case 'none'
        trancorrfn = @(x)x;
    case 'z'
        trancorrfn = @(x)atanh(x);
    otherwise
        error('No matching trancorr parameter');
end

% function used to transform weights
switch I.tranweight
    case 'none'
        tranweightfn = @(x)x;
    case 'sqrt'
        tranweightfn = @(x)sqrt(x);
    case 'sqrtm3'
        tranweightfn = @(x)sqrt(x-3);
    otherwise
        error('No matching tranweight parameter');
end

% similarity function used to compare responses
switch I.simfunc
    case 'corr'
        simfunc = @nanfastcorr;
    case 'nse'
        simfunc = @nan_nse_match_columns;
    case 'mae'
        simfunc = @(a,b)nanmean(-abs(a-b),1);
    otherwise
        error('No matching similarity function')
end

if ~exist(MAT_file, 'file') || I.overwrite
    
    % segment durations
    clear L;
    L.unique_segs = unique(S.segs);
    n_seg_durs = length(L.unique_segs);
    
    % window over which to compute lags
    L.sr = 1/(t(2)-t(1));
    L.lag_t = I.lag_win(1):1/L.sr:I.lag_win(2);
    n_lags = length(L.lag_t);
    
    % number of orders tested
    n_orders = length(unique(S.orders));
    
    % total number of channels to test
    n_channels = length(I.channels);
    
    % number of source stimuli
    n_sourcestim = checkint(S.scramstim_dur/S.sourcestim_dur);
    
    % remove dimensions with only NaNs
    xi = any(any(any(~isnan(D),1),2),3);
    D_noNaN = D(:,:,:,xi);
    n_reps = sum(xi);
    clear xi;
    
    %% Main analysis
    
    L.same_context = zeros(n_lags, n_seg_durs, n_channels);
    L.diff_context = zeros(n_lags, n_seg_durs, n_channels);
    if I.nperms>0
        L.same_context_perm = zeros(n_lags, n_seg_durs, I.nperms, n_channels);
        L.diff_context_perm = zeros(n_lags, n_seg_durs, I.nperms, n_channels);
    end
    L.mean_diff_tc = zeros(n_lags, n_seg_durs, n_channels);
    L.mean_diff_tc_segdur = zeros(n_lags, n_seg_durs, n_seg_durs, n_channels);
    L.n_total_segs = zeros(1,n_seg_durs);
    for q = 1:n_channels
        
        chan = I.channels(q);
        fprintf('\n\nchannel %d\n\n', chan); drawnow;
        
        switch I.boundary
            case 'none'
                n_short_seg_durs = n_seg_durs;
            case {'noleft', 'noright', 'right','left'}
                n_short_seg_durs = n_seg_durs-1;
            case 'nobound'
                n_short_seg_durs = n_seg_durs-2;
            case 'yesbound'
                n_short_seg_durs = n_seg_durs;
            otherwise
                error('No matching boundary constraint');
        end
        for i = 1:n_short_seg_durs
                        
            %% Short segment context
            
            fprintf('\n\nSeg duration: %.0f\n\n', L.unique_segs(i)); drawnow;
            
            % duration and number of segments in the scrambled stimulus
            % and the original source stimuli
            shortseg_dur = L.unique_segs(i)/1000;
            n_shortsegs_per_scramstim = checkint(S.scramstim_dur / shortseg_dur);
            n_shortsegs_per_sourcestim = checkint(S.sourcestim_dur / shortseg_dur);
            
            % onset of each segment in the original source stimulus
            shortseg_onsets_in_sourcestim = (0:n_shortsegs_per_sourcestim-1)*shortseg_dur;
            
            % extract the segment responses for the short segments
            Y_shortseg = nan(n_shortsegs_per_scramstim, n_lags, n_orders, n_reps);
            for l = 1:n_orders
                
                % stimulus corresponding to a particular segment duration and ordering
                xi = S.segs == L.unique_segs(i) & S.orders == l;
                assert(sum(xi)==1);
                
                % data for a given stimulus / channel
                % time x stimulus x channel x repetition
                % -> time x repetition
                X = squeeze_dims(D_noNaN(:,xi,chan,:),[2,3]);
                
                % shortseg_order_indices is n_segs x n_sourcestim matrix
                % and for each segment gives a number that can be used
                % to find the location in the scrambled stimulus by comparing the value
                % with the vector in shortsegs
                shortsegs = S.segorder(xi);
                shortseg_order_indices = reshape(sort(shortsegs.order), n_shortsegs_per_sourcestim, n_sourcestim);
                assert(length(shortsegs.order)==n_shortsegs_per_scramstim);
                clear xi;
                for k = 1:n_shortsegs_per_scramstim % loop through all segments
                    
                    % find location in scrambled stim
                    seg_index_in_scramstim = find(shortseg_order_indices(k)==shortsegs.order);
                    seg_start_time_in_scramstim = (seg_index_in_scramstim-1)*shortseg_dur;
                    
                    % extract response
                    targ_times = L.lag_t' + seg_start_time_in_scramstim;
                    ti = ge_tol(targ_times, t(1)) & le_tol(targ_times, t(end));
                    Y_shortseg(k, ti, l, :) = interp1( t, X, targ_times(ti) );
                    clear seg_index_in_scramstim seg_start_time_in_scramstim;
                end
                clear ti;
                clear X;
            end
            
            %% Longer (or equal) segment responses
            
            % now find responses for corresponding segments in the longer (or equal duration) stimuli
            % boundary constraint controls which segment durations we can use
            switch I.boundary
                case 'none'
                    long_seg_dur_inds = i:n_seg_durs;
                case {'noleft', 'noright', 'right','left'}
                    long_seg_dur_inds = i+1:n_seg_durs;
                case 'nobound'
                    long_seg_dur_inds = i+2:n_seg_durs;
                case 'yesbound'
                    long_seg_dur_inds = i;
                otherwise
                    error('No matching boundary constraint');
            end
            if ~isnan(I.longsegdurs)
                inds_to_select = find(eq_tol(L.unique_segs, I.longsegdurs, 'tol', 1e-6));
                assert(length(inds_to_select)==length(I.longsegdurs));
                long_seg_dur_inds = intersect(long_seg_dur_inds, inds_to_select);
                clear inds_to_select;
            end
                
            %%
            n_long_seg_durs = length(long_seg_dur_inds);
            Y_longsegs = nan(n_shortsegs_per_scramstim, n_lags, n_orders, n_long_seg_durs, n_reps);
            valid_segs = false(n_shortsegs_per_scramstim, n_orders, n_long_seg_durs);
            for j = 1:n_long_seg_durs
                
                fprintf('Responses for seg dur %.0f ms\n', L.unique_segs(long_seg_dur_inds(j))); drawnow;
                
                % duration and number of segments in the scrambled stim
                % and the original source stimuli
                longseg_dur = L.unique_segs(long_seg_dur_inds(j))/1000;
                n_longsegs_per_scramstim = checkint(S.scramstim_dur / longseg_dur);
                n_longsegs_per_sourcestim = checkint(S.sourcestim_dur / longseg_dur);
                
                % onset of each segment in the original source stimulus
                longseg_onsets_in_sourcestim = (0:n_longsegs_per_sourcestim-1)*longseg_dur;
                
                % extract the segment responses for the longer segments
                for l = 1:n_orders
                    
                    % pick out stimulus corresponding to a particular
                    % segment duration and order
                    xi = S.segs == L.unique_segs(long_seg_dur_inds(j)) & S.orders == l;
                    
                    % data for this stimulus/channel
                    % time x stim x channel x repetition
                    % -> time x repetition
                    X = squeeze_dims(D_noNaN(:,xi,chan,:),[2,3]);
                    
                    % longsegs_order_indices is n_segs x n_sourcestim
                    % and for each segment gives a number that can be used
                    % to find the location in the scrambled stim by comparing the value
                    % with the vector in longsegs.order
                    longsegs = S.segorder(xi);
                    longsegs_order_indices = reshape(sort(longsegs.order), n_longsegs_per_sourcestim, n_sourcestim);
                    assert(length(longsegs.order)==n_longsegs_per_scramstim);
                    clear xi;
                    
                    % loop through segments
                    for k = 1:n_shortsegs_per_scramstim
                        
                        % find the start time of the segment in the source stimulus
                        % using this, find which longer segment the shorter segment was part of
                        % and the onset time of the short segment in the longer segment
                        [shortseg_index_in_sourcestim, sourcestim_index] = ind2sub(size(shortseg_order_indices), k);
                        y = shortseg_onsets_in_sourcestim(shortseg_index_in_sourcestim) - longseg_onsets_in_sourcestim;
                        y(y < -1e-4) = inf;
                        [~, which_longseg] = min(y);
                        shortseg_onset_relative_to_longseg = y(which_longseg);
                        clear y;
                        
                        % find if shorter segment is on the boundary of the longer segment
                        on_left_boundary = abs(shortseg_onset_relative_to_longseg-0)<1e-4; % eq_tol(shortseg_onset_relative_to_longseg, 0, 'tol', 1e-4);
                        on_right_boundary = abs(shortseg_onset_relative_to_longseg-(longseg_dur - shortseg_dur))<1e-4; %eq_tol(shortseg_onset_relative_to_longseg, longseg_dur - shortseg_dur, 'tol', 1e-4);
                        on_either_boundary = on_left_boundary || on_right_boundary;
                        on_both_boundary = on_left_boundary && on_right_boundary;
                        switch I.boundary
                            case 'none'
                                valid_segs(k, l, j) = true;
                            case 'noleft'
                                if ~on_left_boundary
                                    valid_segs(k, l, j) = true;
                                else
                                    valid_segs(k, l, j) = false;
                                end
                            case 'noright'
                                if ~on_right_boundary
                                    valid_segs(k, l, j) = true;
                                else
                                    valid_segs(k, l, j) = false;
                                end
                            case 'right'
                                if on_right_boundary
                                    valid_segs(k, l, j) = true;
                                else
                                    valid_segs(k, l, j) = false;
                                end
                            case 'left'
                                if on_left_boundary
                                    valid_segs(k, l, j) = true;
                                else
                                    valid_segs(k, l, j) = false;
                                end
                            case 'nobound'
                                if ~on_either_boundary
                                    valid_segs(k, l, j) = true;
                                else
                                    valid_segs(k, l, j) = false;
                                end
                            case 'yesbound'
                                if on_both_boundary
                                    valid_segs(k, l, j) = true;
                                else
                                    valid_segs(k, l, j) = false;
                                end
                            otherwise
                                error('No matching boundary constraint');
                        end
                        
                        % if a valid segment based on boundary conditions assign
                        if valid_segs(k, l, j)
                            
                            % index of the longer segment int he scrambled stim
                            longseg_index_in_scramstim = find(longsegs_order_indices(which_longseg, sourcestim_index)==longsegs.order);
                            
                            % find corresponding onset time
                            longseg_start_time = (longseg_index_in_scramstim-1)*longseg_dur;
                            
                            % add time to the start of the short segment
                            shortseg_start_time = longseg_start_time + shortseg_onset_relative_to_longseg;
                            
                            % pick out response
                            targ_times = L.lag_t' + shortseg_start_time;
                            ti = targ_times>(t(1)-1e-6) & targ_times < (targ_times+1e-6); % ge_tol(targ_times, t(1)) & le_tol(targ_times, t(end));
                            Y_longsegs(k, ti, l, j, :) = interp1( t, X, targ_times(ti), 'linear' );
                        end
                    end
                    clear ti targ_times;
                    clear X;
                end
            end
            
            % check natcontext does not depend on order and then remove
            X = bsxfun(@minus, valid_segs(:,1,:), valid_segs);
            assert(all(X(:)==0));
            clear X;
            valid_segs = squeeze_dims(valid_segs(:,1,:),2);
            
            % exclude particular central segments
            stim_labels_shortsegs = repmat(1:n_sourcestim, n_shortsegs_per_sourcestim, 1);
            if ~isempty(I.excludesources)
                valid_segs(ismember(stim_labels_shortsegs(:), I.excludesources(:)),:) = false;
            end
            
            %% Weighted correlation
            
            % LEFT OFF HERE
            
            % correlation for same vs. different order
            weight_sum_diff_context = 0;
            weight_sum_same_context = 0;
            weight_sum_diff_context_segdur = zeros(1, n_long_seg_durs);
            for j = 1:n_long_seg_durs
                
                fprintf('Correlations for seg dur %.0f ms\n', L.unique_segs(long_seg_dur_inds(j))); drawnow;
                valid_seg_single_segdur = valid_segs(:,j);
                n_valid_segs = sum(valid_seg_single_segdur);
                weight = tranweightfn(sum(valid_seg_single_segdur)); % standard error
                if q == 1
                    L.n_total_segs(i) = L.n_total_segs(i) + n_valid_segs;
                end
                for k = 1:n_reps
                    for l = (k+1):n_reps
                        for m = 1:n_orders
                            
                            short_rep1 = Y_shortseg(valid_seg_single_segdur, :, m, k);
                            short_rep2 = Y_shortseg(valid_seg_single_segdur, :, m, l);
                            long_rep1 = Y_longsegs(valid_seg_single_segdur, :, m, j, k);
                            long_rep2 = Y_longsegs(valid_seg_single_segdur, :, m, j, l);
                            
                            if I.simoffset
                                offset = (mean(short_rep1(:)) + mean(short_rep2(:)));
                                short_rep1 = short_rep1 + offset;
                                short_rep2 = short_rep2 + offset;
                            end
                            
                            if long_seg_dur_inds(j)==i
                                orders_to_use = setdiff(1:n_orders, m);
                            else
                                orders_to_use = 1:n_orders;
                            end
                            
                            for p = 0:I.nperms
                                
                                if p > 1
                                    short_rep1 = short_rep1(randperm(n_valid_segs),:);
                                    short_rep2 = short_rep2(randperm(n_valid_segs),:);
                                    long_rep1 = long_rep1(randperm(n_valid_segs),:);
                                    long_rep2 = long_rep2(randperm(n_valid_segs),:);
                                end
                                
                                % note: only doing test-retest here for order m, but will hit all orders
                                if I.interleave
                                    X1 = [short_rep1(1:2:end,:); long_rep2(2:2:end,:)];
                                    X2 = [short_rep2(1:2:end,:); long_rep1(2:2:end,:)];
                                    Y1 = [long_rep2(1:2:end,:); short_rep1(2:2:end,:)];
                                    Y2 = [long_rep1(1:2:end,:); short_rep2(2:2:end,:)];
                                    C = trancorrfn(simfunc(X1, X2))/2 + trancorrfn(simfunc(Y1, Y2))/2;
                                    clear X1 X2 Y1 Y2;
                                else
                                    C = (trancorrfn(simfunc(short_rep1, short_rep2))/2 + trancorrfn(simfunc(long_rep1, long_rep2))/2);
                                end
                                if all(~isnan(C))
                                    if p > 1
                                        L.same_context_perm(:,i,p,q) = L.same_context_perm(:,i,p,q) + C' * weight;
                                    else
                                        L.same_context(:,i,q) = L.same_context(:,i,q) + C' * weight;
                                        weight_sum_same_context = weight_sum_same_context + weight;
                                    end
                                end
                                clear C;

                                for n = orders_to_use
                                    
                                    long_rep1 = Y_longsegs(valid_seg_single_segdur, :, n, j, k);
                                    long_rep2 = Y_longsegs(valid_seg_single_segdur, :, n, j, l);
                                    if p > 1
                                        long_rep1 = long_rep1(randperm(n_valid_segs),:);
                                        long_rep2 = long_rep2(randperm(n_valid_segs),:);
                                    end
                                    
                                    if I.interleave
                                        X1 = [short_rep1(1:2:end,:); long_rep2(2:2:end,:)];
                                        X2 = [long_rep2(1:2:end,:); short_rep1(2:2:end,:)];
                                        Y1 = [short_rep2(1:2:end,:); long_rep1(2:2:end,:)];
                                        Y2 = [long_rep1(1:2:end,:); short_rep2(2:2:end,:)];
                                        C = trancorrfn(simfunc(X1, X2))/2 + trancorrfn(simfunc(Y1, Y2))/2;
                                        clear X1 X2 Y1 Y2;
                                    else
                                        C = trancorrfn(simfunc(short_rep1, long_rep2))/2 + trancorrfn(simfunc(short_rep2, long_rep1))/2;
                                    end
                                    if all(~isnan(C))
                                        if p > 1
                                            L.diff_context_perm(:,i,p,q) = L.diff_context(:,i,p,q) + C' * weight;
                                        else
                                            y = mean(Y_longsegs(valid_seg_single_segdur, :, n, j, l) - Y_shortseg(valid_seg_single_segdur, :, m, k),1);
                                            y = squeeze_dims(y,1);
                                            L.mean_diff_tc_segdur(:,i,long_seg_dur_inds(j),q) = L.mean_diff_tc_segdur(:,i,long_seg_dur_inds(j),q) + y * weight;
                                            L.mean_diff_tc(:,i,q) = L.mean_diff_tc(:,i,q) + y * weight;
                                            L.diff_context(:,i,q) = L.diff_context(:,i,q) + C' * weight;
                                            weight_sum_diff_context = weight_sum_diff_context + weight;
                                            weight_sum_diff_context_segdur(j) = weight_sum_diff_context_segdur(j) + weight;
                                        end
                                    end
                                    clear C;
                                end
                            end
                        end
                    end
                end
                L.mean_diff_tc_segdur(:,i,long_seg_dur_inds(j),q) = L.mean_diff_tc_segdur(:,i,long_seg_dur_inds(j),q) / weight_sum_diff_context_segdur(j);
            end
            
            L.same_context(:,i,q) = L.same_context(:,i,q)/weight_sum_same_context;
            L.diff_context(:,i,q) = L.diff_context(:,i,q)/weight_sum_diff_context;
            L.mean_diff_tc(:,i,q) = L.mean_diff_tc(:,i,q)/weight_sum_diff_context;
            if I.nperms>0
                L.same_context_perm(:,i,:,q) = L.same_context_perm(:,i,:,q)/weight_sum_same_context;
                L.diff_context_perm(:,i,:,q) = L.diff_context_perm(:,i,:,q)/weight_sum_diff_context;
            end
            
        end
    end
    
    
    save(MAT_file, 'L');
    
else
    
    load(MAT_file, 'L');
    
end

% save additional parameters
L.channels = I.channels;
L.param_string = param_string;
L.chnames = I.chnames;
L.boundary = I.boundary;
L.figure_directory = [I.figure_directory '/' param_string];
L.output_directory = [I.output_directory '/' param_string];
L.rampwin = S.rampwin;
L.rampdur = S.rampdur;
save(MAT_file, 'L');

%% Plotting

if I.plot_figure
    
    n_seg_durs = length(L.unique_segs);
    
    switch L.boundary
        case 'none'
            n_short_seg_durs = n_seg_durs;
        case {'noleft', 'noright', 'right','left'}
            n_short_seg_durs = n_seg_durs-1;
        case 'nobound'
            n_short_seg_durs = n_seg_durs-2;
        case 'yesbound'
            n_short_seg_durs = n_seg_durs;
        otherwise
            error('No matching boundary constraint');
    end
    
    n_channels = length(L.channels);
    
    % create a figure handle if not supplied
    if isempty(I.figh)
        I.figh = figure;
    end
    for q = 1:n_channels
        
        chan = L.channels(q);
                
        if strcmp(L.chnames{chan}, ['ch' num2str(chan)])
            chname = L.chnames{chan};
        else
            chname = ['ch' num2str(chan) '-' L.chnames{chan}];
        end
        
        %% adaptation plot
        
        if I.plot_adaptation
            clf(I.figh);
            set(I.figh, 'Position', [100, 100, 900, 900]);
            figure(I.figh);
            X = L.mean_diff_tc_segdur(:,:,:,q);
            adapt_range = quantile(X(:), [0.001, 0.999]);
            for k = 1:n_short_seg_durs
                
                subplot(4, 2, k);
                
                % now find responses for corresponding segments in the longer (or equal duration) stimuli
                % boundary constraint controls which segment durations we can use
                switch I.boundary
                    case 'none'
                        long_seg_dur_inds = k:n_seg_durs;
                    case {'noleft', 'noright', 'right','left'}
                        long_seg_dur_inds = k+1:n_seg_durs;
                    case 'nobound'
                        long_seg_dur_inds = k+2:n_seg_durs;
                    case 'yesbound'
                        long_seg_dur_inds = k;
                    otherwise
                        error('No matching boundary constraint');
                end
                n_long_seg_durs = length(long_seg_dur_inds);
                
                h = plot(L.lag_t * 1000, squeeze_dims(L.mean_diff_tc_segdur(:,k,long_seg_dur_inds,q),2), 'LineWidth', 2);
                hold on;
                plot([0 0], adapt_range, 'r--', 'LineWidth', 2);
                plot(L.unique_segs(k)*[1 1], adapt_range, 'r--', 'LineWidth', 2);
                plot(I.plot_win * 1000, [0 0], 'k--', 'LineWidth', 2);
                xlim(I.plot_win * 1000);
                ylim(adapt_range);
                xlabel('Time Lag (ms)');
                ylabel('Mean tc diff');
                title(sprintf('Seg: %.0f ms', L.unique_segs(k)))
                legend_names = cell(1, n_long_seg_durs);
                for l = 1:n_long_seg_durs
                    legend_names{l} = sprintf('%.0f ms', L.unique_segs(long_seg_dur_inds(l)));
                end
                legend(h, legend_names{:}, 'Location', 'Best');
                
            end
            export_fig(mkpdir([L.figure_directory '/cross-context-corr/' chname ...
                '-win' num2str(I.plot_win(1)) '-' num2str(I.plot_win(2)) '-adapt' ...
                '-range' num2str(adapt_range(1), '%.2f') '-' num2str(adapt_range(2), '%.2f') '.pdf']), '-pdf', '-transparent');
        end
        
        %% Lag results
        
        clf(I.figh);
        set(I.figh, 'Position', [100, 100, 900, 900]);
        if isnan(I.plot_range)
            ti = L.lag_t >= I.plot_win(1) & L.lag_t <= I.plot_win(2);
            X = cat(3, L.same_context(ti,:,q), L.diff_context(ti,:,q));
            if any(ismember({'mae'}, I.simfunc))
                corr_range = quantile(X(X~=0), [0.01, 0.99]);
            else
                corr_range = [-1 1] * max(X(:))*1.05;
            end
            clear ti X;
        else
            corr_range = [0 1];
        end
        for k = 1:n_short_seg_durs
            subplot(4, 2, k);
            X = [L.same_context(:,k,q), L.diff_context(:,k,q)];
            plot(L.lag_t * 1000, X, 'LineWidth', 2);
            hold on;
            plot([0 0], corr_range, 'r--', 'LineWidth', 2);
            plot(L.unique_segs(k)*[1 1], corr_range, 'r--', 'LineWidth', 2);
            plot(I.plot_win * 1000, [0 0], 'k--', 'LineWidth', 2);
            xlim(I.plot_win * 1000);
            ylim(corr_range);
            xlabel('Time Lag (ms)');
            ylabel('Z-trans Corr');
            title(sprintf('Seg: %.0f ms', L.unique_segs(k)))
            if k == 1
                legend('Same', 'Diff', 'Diff (Control)');
            end
        end
        export_fig(mkpdir([L.figure_directory '/cross-context-corr/' chname ...
            '-win' num2str(I.plot_win(1)) '-' num2str(I.plot_win(2)) ...
            '-range' num2str(corr_range(1), '%.2f') '-' num2str(corr_range(2), '%.2f') '.pdf']), '-pdf', '-transparent');
    end
end