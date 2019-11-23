function M = modelfit_cross_context_corr_wrapper(...
    L, modelfit_params, output_directory, figure_directory, varargin)

% The dominant computational cost for model fitting is the number of model windows tested not
% the number of channels. So this function combines data across multiple
% lagged correlation analyses to reduce computational cost.
% 
% 2019-11-22: Reshape back to input, plus other mods, Sam NH

dims = size(L);
L = L(:);

% number of lagged correlation structures to combine across
n_groups = length(L);

% default name for each structure
I.groupnames = cell(1, n_groups);
for i = 1:n_groups
    I.groupnames{i} = ['group' num2str(i)];
end

% overwrite with optional inputs
I = parse_optInputs_keyvalue(varargin, I);

% combine
L_all = L{1};
L_all.chnames = strcat(I.groupnames{1}, '-', L{1}.chnames);
for i = 2:n_groups
    
    % combine correlation, combine over channel dimension
    L_all.same_context = cat(3, L_all.same_context, L{i}.same_context);
    L_all.diff_context = cat(3, L_all.diff_context, L{i}.diff_context);
    
    % misc adaptation parameters, combine over channel dimension
    L_all.mean_diff_tc = cat(3, L_all.mean_diff_tc, L{i}.mean_diff_tc);
    L_all.mean_diff_tc_segdur = cat(4, L_all.mean_diff_tc_segdur, L{i}.mean_diff_tc_segdur);
    
    % combine across channels and group names, add group name to the channel name
    L_all.chnames = [L_all.chnames, strcat(I.groupnames{i}, '-', L{i}.chnames)];
    L_all.channels = cat(2, L_all.channels, L{i}.channels);
    
    % check to make sure they have the same parameters
    f_to_check = {'lag_t', 'sr', 'unique_segs', 'rampwin', 'rampdur'};
    for j = 1:length(f_to_check)
        assert(isequal(L_all.(f_to_check{j}), L{i}.(f_to_check{j})));
    end
end

% run
L_all.figure_directory = figure_directory;
L_all.output_directory = output_directory;
M_all = modelfit_cross_context_corr(L_all, modelfit_params{:});

% split back up
M = cell(dims);
for i = 1:length(L)
    index = 0;
    n_channels = length(L{i}.channels);
    chan = (1:n_channels) + index;
    M{i} = M_all;
    M{i}.Y = M_all.Y(:,:,chan);
    M{i}.loss = M_all.loss(:,:,:,chan);
    M{i}.Ybest = M_all.Ybest(:,:,chan);
    M{i}.best_intper_sec = M_all.best_intper_sec(chan,:);
    M{i}.best_delay_sec = M_all.best_delay_sec(chan,:);
    M{i}.best_shape = M_all.best_shape(chan,:);
    M{i}.best_loss = M_all.best_loss(chan,:);
    M{i}.channels = M_all.channels(chan);
    if isfield(M_all, 'logP_gaussfit')
        M{i}.logP_gaussfit = M_all.logP_gaussfit(chan);
    end
    if isfield(M_all, 'logP_counts')
        M{i}.logP_counts = M_all.logP_counts(chan);
    end
end
