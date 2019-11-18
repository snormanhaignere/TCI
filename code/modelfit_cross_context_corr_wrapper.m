function M_all = modelfit_cross_context_corr_wrapper(...
    L, groupnames, figure_directory, output_directory, varargin)

% The dominant computational cost for model fitting is the number of model windows tested not
% the number of channels. So this function combines data across multiple
% lagged correlation analyses to reduce computational cost.

% number of lagged correlation analyses
n_groups = length(L);
assert(n_groups == length(groupnames));

% combine
L_all = L{1};
L_all.chnames = strcat(groupnames{1}, L{1}.chnames);
for i = 2:n_groups
    L_all.same_context = cat(3, L_all.same_context, L{i}.same_context);
    L_all.diff_context = cat(3, L_all.diff_context, L{i}.diff_context);
    L_all.mean_diff_tc = cat(3, L_all.mean_diff_tc, L{i}.mean_diff_tc);
    L_all.mean_diff_tc_segdur = cat(4, L_all.mean_diff_tc_segdur, L{i}.mean_diff_tc_segdur);
    L_all.chnames = [L_all.chnames, strcat(groupnames{i}, L{i}.chnames)];
    L_all.channels = cat(2, L_all.channels, (1:length(L{i}.channels)) + length(L_all.channels));
    
    % check to make sure they have the same parameters
    f_to_check = {'lag_t', 'sr', 'unique_segs', 'rampwin', 'rampdur'};
    for j = 1:length(f_to_check)
        assert(isequal(L_all.(f_to_check{j}), L{i}.(f_to_check{j})));
    end
end

% run
L_all.figure_directory = figure_directory;
L_all.output_directory = output_directory;
M_all = modelfit_cross_context_corr(L_all, varargin{:});

    

