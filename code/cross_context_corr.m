function [L, MAT_file] = cross_context_corr(D, t, S, varargin)

% Calculates cross and same context correlation for different lags
% relative to segment onset.
%
% D: time x stim x channel x repetition
%
% t: vector of time stamps; length(t) must equal size(D,1)
%
% S: structure with stimulsus information
%
% See comments below for optional inputs below (fields of I)

%% Optional parameters and defaults

clear I;

% window used to compute lags
I.lag_win = 2*max(S.segs(:))/1000*[0,1];

% window used for plotting
I.plot_win = I.lag_win;

% boundary constraint
% 'any': all segment comparisons allowed
% 'samedur': segments must be of equal duration (i.e. no embedded segs)
% 'diffdur': segments must of different durations (i.e. only embedded segs)
% 'noleft': shorter seg cannot be on left boundary of longer seg
% 'noright': shorter seg cannot be on right boundary of longer seg
% 'noleftright': cannot be on either boundary
% 'right': shorter seg must be on right of boundary
% 'left': shorter seg must be on left of boundary
I.boundary = 'any';

% metric used to compare responses across contexts:
% 'corr': pearson correlation
% 'nse': normalized squared error
% 'mae': minimum absolute error
I.simfunc = 'corr';

% channels for which to perform analysis (default: all)
I.channels = 1:size(D,3);

% whether to interleave shorter and longer segments
I.interleave_diffdur = true;

% whether to interleave segments of different orders but the same duration
I.interleave_samedur = false;

% can optionally force the longer of the two segment durations being
% compared to have a fixed set of values
I.longsegdurs = NaN;

% whether to compare different contexts from the same repetition
I.samerep = false;

% whether to average across odd and even runs
% alternative is to compare all pairs of runs (the default)
% can optionally force there to be an equal number of odd and even runs by
% discarding extra odd runs
I.oddeven = false;
I.oddevenmatch = false;

% can exclude segments that come from particular sources
% the sources are specified in S.sources
% specify the indices of the sources to exclude
I.excludesources = [];

% whether to overwrite or debug
I.keyboard = false;
I.overwrite = false;

% function to transfrom correlation values:
% 'none': no transformation (default)
% 'z': fisher z-transform
I.trancorr = 'none';

% function used to transform weights:
% 'none': no transformation
% 'sqrt': square root
% 'sqrtm3': square root minus three (based on standard error of z values)
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

% number of permutations
I.nperms = 0;

% number of bootstrapped samples
% and type of bootstrapping to perform
% ** Still needs to be coded **
I.nbstraps = 0;
I.bstraptype = 'sources';

% number of splits of segments
I.nsplits = 1;
I.splitfac = 1;
I.splitrand = false;
I.splitbysource = false;

% whether to plot additional
% samples computed via permutation test
% or bootstrapping
I.plot_extra_smps = false;

% random seed, only relevant for random processes
% (e.g. permutation testing or bootstrapping)
I.randseed = 1;

% whether to plot figures
I.plot_figure = true;

% range of values to plot
% if NaN, set to central 98% of plotted values
I.plot_range = NaN;

% figure handle
I.figh = matlab.ui.Figure.empty;

% if false, no analysis is performed, but MAT file is returned
I.run = true;

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

% string with key parameter values
% always include boundary and lag_win
% for other parameters include them if
% they changed from their default values
always_include = {'boundary','lag_win'};
always_exclude = {...
    'plot_win', 'plot_figure', 'plot_range', ...
    'keyboard', 'overwrite', 'output_directory', 'figure_directory', ...
    'chnames', 'run', 'figh', 'plot_extra_smps'};
param_string = optInputs_to_string(I, C_value, always_include, always_exclude);

% MAT file to save results to
MAT_file = mkpdir([I.output_directory '/' param_string '/corr_data.mat']);

% can stop and just return the MAT_file
% (useful to check if it exists)
if ~I.run
    clear L;
    L.param_string = param_string;
    L.figure_directory = [I.figure_directory '/' param_string];
    L.output_directory = [I.output_directory '/' param_string];
    return;
end


%% Do the analysis!

if ~exist(MAT_file, 'file') || I.overwrite
    
    ResetRandStream2(I.randseed);
    
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
    
    % optionally average odd and even reps
    if I.oddeven
        odd_reps = 1:2:n_reps;
        even_reps = 2:2:n_reps;
        if I.oddevenmatch
            odd_reps = odd_reps(1:length(even_reps));
        end
        D_noNaN = cat(4, ...
            nanmean(D_noNaN(:,:,:,odd_reps),4), ...
            nanmean(D_noNaN(:,:,:,even_reps),4));
        n_reps = 2;
        clear odd_reps even_reps;
    end
    
    % number of additional bootstrapped or permuted samples (not both)
    assert(sum(I.nperms>0 && I.nbstraps>0 && I.nsplits) < 2)
    n_smps = max([I.nperms, I.nbstraps, I.nsplits]);
    
    %% Short functions used later on
    
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
    
    % whether to make any same duration comparisons
    make_samedur_comparisons = n_orders > 1 && ~any(ismember({'diffdur', 'noleftright', 'noleft', 'noright'}, I.boundary));
    
    % whether to make any different duration comparisons
    make_diffdur_comparisons = ~any(ismember({'samedur'}, I.boundary));
    
    %% Pairs of orders for same duration comparisons
    
    % for same-duration comparisons, must always compare different orders
    if make_samedur_comparisons
        samedur_order_pairs = [];
        for o1 = 1:n_orders
            for o2 = o1+1:n_orders
                new_pair = [o1, o2]';
                samedur_order_pairs = cat(2, samedur_order_pairs, new_pair);
            end
        end
        clear o1 o2 new_pair;
    end
    
    %% Pairs of orders for different duration comparisons
    
    % since the contexts are always different, we can consider
    % all possible pairs of orders for those different durations
    % i.e. o1-o1 is fine, o1-o2 is different from o2-o1
    if make_diffdur_comparisons
        diffdur_order_pairs = [];
        for o1 = 1:n_orders
            for o2 = 1:n_orders
                new_pair = [o1, o2]';
                diffdur_order_pairs = cat(2, diffdur_order_pairs, new_pair);
            end
        end
        clear o1 o2 new_pair;
    end
    
    %% Pairs of repetitions for same vs. different context comparisons
    
    % pairs of repetitions to use for diff context comparisons
    % when comparing different contexts/order
    % we need to consider all possible pairs of repetitions
    % since order1-rep1 vs. order2-rep2 and order1-rep2 vs. order2-rep1
    % are both valid
    % we can optionally also include pairs from the same
    % repetition, though it is arguably more conservative to not do this
    diffcontext_rep_pairs = [];
    for r1 = 1:n_reps
        if I.samerep
            second_reps = 1:n_reps;
        else
            second_reps = [1:r1-1,r1+1:n_reps];
        end
        for r2 = second_reps
            new_pair = [r1, r2]';
            diffcontext_rep_pairs = cat(2, diffcontext_rep_pairs, new_pair);
        end
    end
    clear r1 r2 second_reps new_pair;
    
    % pairs of repetitions to use for same context comparisons
    % here we need to consider unique rep pairs obviously
    % since we have exactly the same segment/context
    samecontext_rep_pairs = [];
    for r1 = 1:n_reps
        for r2 = r1+1:n_reps
            new_pair = [r1, r2]';
            samecontext_rep_pairs = cat(2, samecontext_rep_pairs, new_pair);
        end
    end
    clear r1 r2 new_pair;
        
    %% Segment information
    
    % Key variables:
    %
    % seg_start_time_in_scramstim
    % {segment duration} x seg x order
    % indicates the onset time for each segment in a scrambled stimulus
    %
    % embedded_seg_start_time_in_scramstim:
    % {segment duration} x seg x order x long seg duration
    % indicates the onset time of a segment embedded in a longer segment
    % segment in a scrambled stimulus
    %
    % samedur_valid_segs
    % {segment duration} x seg
    % indicates whether a segment should be used same-duration comparisons
    %
    % diffdur_valid_segs
    % {segment duration} x seg x long seg duration
    % indicates whether a segment should be used different-duration comparisons
    
    n_segs_per_scramstim = nan(1, n_seg_durs);
    seg_start_time_in_scramstim = cell(1, n_seg_durs);
    embedded_seg_start_time_in_scramstim = cell(1, n_seg_durs);
    samedur_valid_segs = cell(1, n_seg_durs);
    diffdur_valid_segs = cell(1, n_seg_durs);
    longer_seg_dur_inds = cell(1, n_seg_durs);
    n_longer_seg_durs = nan(1, n_seg_durs);
    source_labels = cell(1, n_seg_durs);
    for i = 1:n_seg_durs
        
        % duration and number of segments in the scrambled stimulus
        % and the original source stimuli
        seg_dur = L.unique_segs(i)/1000;
        n_segs_per_scramstim(i) = checkint(S.scramstim_dur / seg_dur);
        n_segs_per_sourcestim = checkint(S.sourcestim_dur / seg_dur);
        
        % shortseg_order_indices is n_segs_per_source_stim x n_sourcestim matrix
        % Each element contains an index that corresponds to a single segment.
        % The location of this index S.segorder.order indicates where this
        % segment occured in the scrambled stimulus.
        xi = find(S.segs == L.unique_segs(i),1);
        seg_order_indices = reshape(sort(S.segorder(xi).order), n_segs_per_sourcestim, n_sourcestim);
        assert(numel(seg_order_indices)==n_segs_per_scramstim(i));
        clear xi;
        
        % for each segment a label indicating the source it came from
        source_labels{i} = repmat(1:n_sourcestim, n_segs_per_sourcestim, 1);
        assert(numel(source_labels{i})==n_segs_per_scramstim(i));
        
        % onset of each segment in the original source stimulus
        seg_onsets_in_sourcestim = (0:n_segs_per_sourcestim-1)*seg_dur;
        seg_start_time_in_scramstim{i} = nan(n_segs_per_scramstim(i), n_orders);
        for l = 1:n_orders
            xi = S.segs == L.unique_segs(i) & S.orders == l;
            assert(sum(xi)==1);
            for k = 1:n_segs_per_scramstim(i) % loop through all segments
                seg_index_in_scramstim = find(seg_order_indices(k)==S.segorder(xi).order);
                seg_start_time_in_scramstim{i}(k,l) = (seg_index_in_scramstim-1)*seg_dur;
            end
            clear xi;
        end
        
        % determine which segments are valid for same-duration comparisons
        % in this case the only reason a segment would be invalid is if 
        % its source was excluded
        if make_samedur_comparisons
            samedur_valid_segs{i} = ~ismember(source_labels{i}(:), I.excludesources(:));
        end
        
        %% Now find embedded segments
        
        % only need to do this if making different duration comparisons
        if make_diffdur_comparisons 
            
            % segment durations to compare with
            longer_seg_dur_inds{i} = i+1:n_seg_durs;
            
            % select particular long segment durations to use for reference
            if ~isnan(I.longsegdurs)
                longsegs_to_use = find(eq_tol(L.unique_segs, I.longsegdurs, 'tol', 1e-6));
                assert(length(longsegs_to_use)==length(I.longsegdurs));
                longer_seg_dur_inds{i} = intersect(longer_seg_dur_inds{i}, longsegs_to_use);
                clear longsegs_to_use;
            end
            n_longer_seg_durs(i) = length(longer_seg_dur_inds{i});
            
            % double check they are sorted
            assert(all(longer_seg_dur_inds{i} == sort(longer_seg_dur_inds{i})));
            
            % if there are any longer segs
            % let's figure out their onset time
            % and whether they are valid given boundary constraints
            if n_longer_seg_durs(i)>0
                
                embedded_seg_start_time_in_scramstim{i} = nan(n_segs_per_scramstim(i), n_orders, n_longer_seg_durs(i));
                diffdur_valid_segs{i} = false(n_segs_per_scramstim(i), n_longer_seg_durs(i));
                for j = 1:n_longer_seg_durs(i)
                    
                    % duration and number of segments in the scrambled stim
                    % and the original source stimuli
                    longer_seg_dur = L.unique_segs(longer_seg_dur_inds{i}(j))/1000;
                    n_longer_segs_per_scramstim = checkint(S.scramstim_dur / longer_seg_dur);
                    n_longer_segs_per_sourcestim = checkint(S.sourcestim_dur / longer_seg_dur);
                    
                    % onset of each segment in the original source stimulus
                    longer_seg_onsets_in_sourcestim = (0:n_longer_segs_per_sourcestim-1)*longer_seg_dur;
                    
                    % longer_seg_order_indices is n_segs_per_source_stim x n_sourcestim matrix
                    % Each element contains an index that corresponds to a single segment.
                    % The location of this index S.segorder.order indicates where this
                    % segment occured in the scrambled stimulus.
                    xi = find(S.segs == L.unique_segs(longer_seg_dur_inds{i}(j)),1);
                    longer_seg_order_indices = reshape(sort(S.segorder(xi).order), n_longer_segs_per_sourcestim, n_sourcestim);
                    assert(numel(longer_seg_order_indices)==n_longer_segs_per_scramstim);
                    clear xi;
                    
                    % loop through segments
                    for k = 1:n_segs_per_scramstim(i)
                        
                        % if seg should be excluded, we can skip the rest
                        if any(source_labels{i}(k)==I.excludesources)
                            
                            diffdur_valid_segs{i}(k, j) = false;
                            
                        else
                            
                            % find the start time of the shorter segment in the source stimulus
                            % using this, find which longer segment the shorter segment was part of
                            % and the onset time of the short segment in the longer segment
                            [seg_index_in_sourcestim, sourcestim_index] = ind2sub(size(seg_order_indices), k);
                            y = seg_onsets_in_sourcestim(seg_index_in_sourcestim) - longer_seg_onsets_in_sourcestim;
                            y(y < -1e-4) = inf;
                            [~, which_longseg] = min(y);
                            seg_onset_relative_to_longer_seg = y(which_longseg);
                            clear y;
                            
                            % find if shorter segment is on the boundary of the longer segment
                            on_left_boundary = abs(seg_onset_relative_to_longer_seg-0)<1e-6; % eq_tol(shortseg_onset_relative_to_longseg, 0, 'tol', 1e-4);
                            on_right_boundary = abs(seg_onset_relative_to_longer_seg-(longer_seg_dur - seg_dur))<1e-6; %eq_tol(shortseg_onset_relative_to_longseg, longseg_dur - shortseg_dur, 'tol', 1e-4);
                            on_either_boundary = on_left_boundary || on_right_boundary;
                            switch I.boundary
                                case {'any', 'diffdur'}
                                    diffdur_valid_segs{i}(k, j) = true;
                                case 'noleft'
                                    if ~on_left_boundary
                                        diffdur_valid_segs{i}(k, j) = true;
                                    else
                                        diffdur_valid_segs{i}(k, j) = false;
                                    end
                                case 'noright'
                                    if ~on_right_boundary
                                        diffdur_valid_segs{i}(k, j) = true;
                                    else
                                        diffdur_valid_segs{i}(k, j) = false;
                                    end
                                case 'right'
                                    if on_right_boundary
                                        diffdur_valid_segs{i}(k, j) = true;
                                    else
                                        diffdur_valid_segs{i}(k, j) = false;
                                    end
                                case 'left'
                                    if on_left_boundary
                                        diffdur_valid_segs{i}(k, j) = true;
                                    else
                                        diffdur_valid_segs{i}(k, j) = false;
                                    end
                                case 'noleftright'
                                    if ~on_either_boundary
                                        diffdur_valid_segs{i}(k, j) = true;
                                    else
                                        diffdur_valid_segs{i}(k, j) = false;
                                    end
                                otherwise
                                    error('No matching boundary constraint');
                            end
                            
                            % if a valid segment based on boundary conditions
                            % then figure out the onset time relative to the scrambled stimulus
                            if diffdur_valid_segs{i}(k, j)
                                
                                % extract the segment responses for the longer segments
                                for l = 1:n_orders
                                    
                                    % index of the longer segment int he scrambled stim
                                    xi = S.segs == L.unique_segs(longer_seg_dur_inds{i}(j)) & S.orders == l;
                                    assert(sum(xi)==1);
                                    longer_seg_index_in_scramstim = find(longer_seg_order_indices(which_longseg, sourcestim_index)==S.segorder(xi).order);
                                    assert(length(longer_seg_index_in_scramstim)==1);
                                    clear xi;
                                    
                                    % find corresponding onset time
                                    longer_seg_start_time = (longer_seg_index_in_scramstim-1)*longer_seg_dur;
                                    
                                    % add time to the start of the short segment
                                    embedded_seg_start_time_in_scramstim{i}(k,l,j) = longer_seg_start_time + seg_onset_relative_to_longer_seg;
                                    clear longer_seg_start_time;
                                end
                            end
                        end
                    end
                end
            end
        end
                
        %% double check the excluded sources are absent
        
        if ~isempty(I.excludesources)
            assert(isempty(intersect(source_labels{i}(samedur_valid_segs{i}), I.excludesources)));
            for j = 1:n_longer_seg_durs(i)
                assert(isempty(intersect(source_labels{i}(diffdur_valid_segs{i}(:,j)), I.excludesources)));
            end
        end
        
    end
    
    %% Now actually do the analysis using the above info
    
    fprintf('Correlation computation\n');
    L.same_context_twogroups = zeros(n_lags, n_seg_durs, n_channels, n_smps+1, 2);
    L.same_context = zeros(n_lags, n_seg_durs, n_channels, n_smps+1);
    L.same_context_err = zeros(n_lags, n_seg_durs, n_channels, n_smps+1);
    L.diff_context = zeros(n_lags, n_seg_durs, n_channels, n_smps+1);
    L.n_total_segs = zeros(n_seg_durs,n_smps+1);
    for q = 1:n_channels % analysis is done separately for every channel to save memory
        
        chan = I.channels(q);
        fprintf('\n\nchannel %d\n', chan); drawnow;
        
        for i = 1:n_seg_durs
            
            %% Get the segments
            
            % skip this segment duration if there are no valid segment comparisons
            if all(~diffdur_valid_segs{i}(:)) && all(~samedur_valid_segs{i}(:))
                continue;
            end
            
            fprintf('Seg duration: %.0f\n', L.unique_segs(i)); drawnow;
            
            % non-embedded segments
            Y_seg = nan(n_segs_per_scramstim(i), n_lags, n_orders, n_reps);
            for l = 1:n_orders
                xi = S.segs == L.unique_segs(i) & S.orders == l;
                assert(sum(xi)==1);
                X = squeeze_dims(D_noNaN(:,xi,chan,:),[2,3]);
                for k = 1:n_segs_per_scramstim(i) % loop through all segments
                    if ((make_samedur_comparisons && samedur_valid_segs{i}(k)) || ...
                        (make_diffdur_comparisons && any(diffdur_valid_segs{i}(k, :))))
                        targ_times = L.lag_t' + seg_start_time_in_scramstim{i}(k,l);
                        ti = targ_times > (t(1)-1e-6) & targ_times < (t(end)+1e-6);
                        Y_seg(k, ti, l, :) = interp1( t, X, targ_times(ti) );
                    end
                end
                clear ti targ_times;
            end
            clear X;
            
            % embedded segments, only needed for different duration comparisons
            if make_diffdur_comparisons
                Y_embed_seg = nan(n_segs_per_scramstim(i), n_lags, n_orders, n_reps, n_longer_seg_durs(i));
                for j = 1:n_longer_seg_durs(i)
                    for l = 1:n_orders
                        xi = S.segs == L.unique_segs(longer_seg_dur_inds{i}(j)) & S.orders == l;
                        assert(sum(xi)==1);
                        X = squeeze_dims(D_noNaN(:,xi,chan,:),[2,3]);
                        for k = 1:n_segs_per_scramstim(i)
                            if diffdur_valid_segs{i}(k, j)
                                targ_times = L.lag_t' + embedded_seg_start_time_in_scramstim{i}(k,l,j);
                                ti = targ_times>(t(1)-1e-6) & targ_times < (targ_times+1e-6);
                                Y_embed_seg(k, ti, l, :, j) = interp1( t, X, targ_times(ti), 'linear' );
                            end
                        end
                        clear X ti targ_times;
                    end
                end
            end
            
            %% Split index if we're doing splits
            
            if I.nsplits > 1
                N = n_segs_per_scramstim(i);
                assert(N==length(samedur_valid_segs{i}))
                if I.splitbysource
                    chunk_index = source_labels{i}(:)-1;
                else
                    chunkforsplit = N/(I.nsplits*I.splitfac);
                    chunk_index = floor((0:N-1)/(chunkforsplit));
                    clear chunkforsplit;
                end
                if I.splitrand
                    xi = randperm(chunk_index(end)+1);
                    chunk_index = xi(chunk_index+1)-1;
                    clear xi;
                end
                split_index = mod(chunk_index, I.nsplits)+1;
                clear N chunk_index;
            end
            
            %% Correlations
            
            for b = 1:n_smps+1
                
                same_context_weight_twogroups = [0,0];
                diff_context_weight = 0;
                
                %% Segs for the same duration
                
                if make_samedur_comparisons
                    
                    samedur_seg_inds = find(samedur_valid_segs{i});
                    
                    % optionally bootstrap segments
                    if b > 1 && I.nbstraps > 0
                        error('Need to implement bootstrapping');
                    end
                    
                    % select segments for different splits
                    if b > 1 && I.nsplits > 1
                        if ~isempty(samedur_seg_inds)
                            samedur_seg_inds = intersect(samedur_seg_inds, find(split_index==(b-1)));
                            assert(~isempty(samedur_seg_inds));
                        end
                    end
                    
                    % add to the count
                    if q == 1
                        L.n_total_segs(i,b) = L.n_total_segs(i,b) + length(samedur_seg_inds);
                    end
                    
                    
                    % optionally permute segment orders
                    if b > 1 && I.nperms > 0
                        xi = randperm(length(samedur_seg_inds));
                        samedur_seg_pairs = [samedur_seg_inds, samedur_seg_inds(xi)];
                        clear xi;
                    else
                        samedur_seg_pairs = [samedur_seg_inds, samedur_seg_inds];
                    end
                end
                
                %% Segs for different duration
                
                if make_diffdur_comparisons
                    
                    diffdur_seg_pairs = cell(1, n_longer_seg_durs(i));
                    for j = 1:n_longer_seg_durs(i)
                        
                        diffdur_seg_inds = find(diffdur_valid_segs{i}(:,j));
                        
                        % optionally bootstrap segments
                        if b > 1 && I.nbstraps>0
                            error('Need to implement bootstrapping');
                        end
                        
                        % select segments for different splits
                        if b > 1 && I.nsplits > 1
                            if ~isempty(diffdur_seg_inds)
                                diffdur_seg_inds = intersect(diffdur_seg_inds, find(split_index==(b-1)));
                                assert(~isempty(diffdur_seg_inds));
                            end
                        end
                        
                        % add to the count
                        L.n_total_segs(i,b) = L.n_total_segs(i,b) + length(diffdur_seg_inds);
                        
                        % optionally permute segment orders
                        if b > 1 && I.nperms>0
                            xi = randperm(length(diffdur_seg_inds));
                            diffdur_seg_pairs{j} = [diffdur_seg_inds, diffdur_seg_inds(xi)];
                            clear xi;
                        else
                            diffdur_seg_pairs{j} = [diffdur_seg_inds, diffdur_seg_inds];
                        end
                    end
                end
                
                %% Compare same duration segments
                
                if make_samedur_comparisons
                    
                    % weight by number of segments in the comparison
                    weight = tranweightfn(size(samedur_seg_pairs,1));
                    X = cell(1,2);
                    for k = 1:size(samedur_order_pairs,2)
                        p_orders = samedur_order_pairs(:,k);
                        assert(p_orders(1)~=p_orders(2));
                        
                        % select segs to be compared
                        X{1} = Y_seg(samedur_seg_pairs(:,1), :, p_orders(1), :);
                        X{2} = Y_seg(samedur_seg_pairs(:,2), :, p_orders(2), :);
                        assert(all(size(X{1})==size(X{2})));
                        
                        % optionally interleave values
                        if I.interleave_samedur
                            [X{1}, X{2}] = interleave_oddeven(X{1}, X{2});
                        end
                        
                        % correlate pairs of repetitions
                        for l = 1:size(diffcontext_rep_pairs,2)
                            p_reps = diffcontext_rep_pairs(:,l);
                            C = trancorrfn(simfunc(X{1}(:,:,1,p_reps(1)), X{2}(:,:,1,p_reps(2))));
                            if all(~isnan(C))
                                L.diff_context(:,i,q,b) = L.diff_context(:,i,q,b) + C' * weight;
                                diff_context_weight = diff_context_weight + weight;
                            end
                        end
                        
                        % reliability of each element
                        for l = 1:size(samecontext_rep_pairs,2)
                            p_reps = samecontext_rep_pairs(:,l);
                            for m = 1:2
                                C = trancorrfn(simfunc(X{m}(:,:,1,p_reps(1)), X{m}(:,:,1,p_reps(2))));
                                if all(~isnan(C))
                                    L.same_context_twogroups(:,i,q,b,m) = L.same_context_twogroups(:,i,q,b,m) + C' * weight;
                                    same_context_weight_twogroups(m) = same_context_weight_twogroups(m) + weight;
                                end
                            end
                        end
                    end
                    clear X;
                end
                
                %% Different duration
                
                if make_diffdur_comparisons
                    X = cell(1,2);
                    for j = 1:n_longer_seg_durs(i)
                        weight = tranweightfn(size(diffdur_seg_pairs{j},1)); % standard error
                        for k = 1:size(diffdur_order_pairs,2)
                            p_orders = diffdur_order_pairs(:,k);
                            
                            % select segs to be compared
                            X{1} = Y_seg(diffdur_seg_pairs{j}(:,1), :, p_orders(1), :);
                            X{2} = Y_embed_seg(diffdur_seg_pairs{j}(:,1), :, p_orders(2), :, j);
                            if ~(all(size(X{1})==size(X{2})))
                                keyboard;
                            end
                            
                            % optionally interleave values
                            if I.interleave_diffdur
                                [X{1}, X{2}] = interleave_oddeven(X{1}, X{2});
                            end
                            
                            % correlate pairs of repetitions
                            for l = 1:size(diffcontext_rep_pairs,2)
                                p_reps = diffcontext_rep_pairs(:,l);
                                C = trancorrfn(simfunc(X{1}(:,:,1,p_reps(1)), X{2}(:,:,1,p_reps(2))));
                                if all(~isnan(C))
                                    L.diff_context(:,i,q,b) = L.diff_context(:,i,q,b) + C' * weight;
                                    diff_context_weight = diff_context_weight + weight;
                                end
                            end
                            
                            % reliability of each element
                            for l = 1:size(samecontext_rep_pairs,2)
                                p_reps = samecontext_rep_pairs(:,l);
                                for m = 1:2
                                    C = trancorrfn(simfunc(X{m}(:,:,1,p_reps(1)), X{m}(:,:,1,p_reps(2))));
                                    if all(~isnan(C))
                                        L.same_context_twogroups(:,i,q,b,m) = L.same_context_twogroups(:,i,q,b,m) + C' * weight;
                                        same_context_weight_twogroups(m) = same_context_weight_twogroups(m) + weight;
                                    end
                                end
                            end
                        end
                    end
                    clear X;
                end
                
                %% Divide by weights, combine groups
                
                % divide by weights
                L.same_context_twogroups(:,i,q,b,1) = L.same_context_twogroups(:,i,q,b,1)/same_context_weight_twogroups(1);
                L.same_context_twogroups(:,i,q,b,2) = L.same_context_twogroups(:,i,q,b,2)/same_context_weight_twogroups(2);
                L.diff_context(:,i,q,b) = L.diff_context(:,i,q,b)/diff_context_weight;
                
                % average and substract the two groups
                L.same_context(:,i,q,b) = L.same_context_twogroups(:,i,q,b,1)/2 + L.same_context_twogroups(:,i,q,b,2)/2;
                L.same_context_err(:,i,q,b) = (L.same_context_twogroups(:,i,q,b,1)/2 - L.same_context_twogroups(:,i,q,b,2)/2).^2;
                
            end
        end
    end
    
    % additional parameters
    L.channels = I.channels;
    L.param_string = param_string;
    L.chnames = I.chnames(I.channels);
    L.boundary = I.boundary;
    L.figure_directory = [I.figure_directory '/' param_string];
    L.output_directory = [I.output_directory '/' param_string];
    L.rampwin = S.rampwin;
    L.rampdur = S.rampdur;
    
    % bootstrapped error
    if I.nbstraps>0
        X = L.same_context(:,:,:,2:end);
        Y = bsxfun(@minus, X, mean(X,4)).^2;
        L.same_context_bstrap_err = mean(Y(:,:,:,:),4);
        clear X Y;
    end
    
    save(MAT_file, 'L');
    
else
    
    load(MAT_file, 'L');
    
end

%% Plotting

if I.plot_figure
    
    n_seg_durs = length(L.unique_segs);
    n_channels = length(L.channels);
    n_smps = size(L.same_context,4);
    
    % create a figure handle if not supplied
    if isempty(I.figh)
        I.figh = figure;
    end
    for q = 1:n_channels
        
        % channel to plot
        chan = L.channels(q);
        
        % channel name for saving
        if strcmp(L.chnames{q}, ['ch' num2str(chan)])
            chname = L.chnames{q};
        else
            chname = ['ch' num2str(chan) '-' L.chnames{q}];
        end
        
        % samples to plot
        if I.plot_extra_smps
            smps_to_plot = 1:n_smps;
        else
            smps_to_plot = 1;
        end
        for b = smps_to_plot
            
            % cross and same context correlations
            diff_context = L.diff_context(:,:,:,b);
            same_context = L.same_context(:,:,:,b);
            
            % plot
            clf(I.figh);
            set(I.figh, 'Position', [100, 100, 900, 900]);
            if isnan(I.plot_range)
                ti = L.lag_t >= I.plot_win(1) & L.lag_t <= I.plot_win(2);
                X = cat(3, same_context(ti,:,q), diff_context(ti,:,q));
                if any(ismember({'mae'}, I.simfunc))
                    corr_range = quantile(X(X~=0), [0.01, 0.99]);
                else
                    corr_range = [-1 1] * max(X(:))*1.05;
                end
                clear ti X;
            else
                corr_range = I.plot_range;
            end
            for k = 1:n_seg_durs
                subplot(4, 2, k);
                X = [same_context(:,k,q), diff_context(:,k,q)];
                plot(L.lag_t * 1000, X, 'LineWidth', 2);
                hold on;
                plot([0 0], corr_range, 'r--', 'LineWidth', 2);
                plot(L.unique_segs(k)*[1 1], corr_range, 'r--', 'LineWidth', 2);
                plot(I.plot_win * 1000, [0 0], 'k--', 'LineWidth', 2);
                xlim(I.plot_win * 1000);
                ylim(corr_range);
                xlabel('Time Lag (ms)');
                ylabel('Similarity');
                title(sprintf('Seg: %.0f ms', L.unique_segs(k)))
                if k == 1
                    legend('Same', 'Cross');
                end
            end
            if b > 0
                smpstr = ['-smp' num2str(b)];
            else
                smpstr = '';
            end
            fname = [L.figure_directory '/cross-context-corr/' chname ...
                '-win' num2str(I.plot_win(1)) '-' num2str(I.plot_win(2)) ...
                smpstr '-range' num2str(corr_range(1), '%.2f') '-' num2str(corr_range(2), '%.2f')];

            export_fig(mkpdir([fname '.pdf']), '-pdf', '-transparent');
            export_fig(mkpdir([fname '.png']), '-png', '-transparent', '-r150');
        end
    end
end