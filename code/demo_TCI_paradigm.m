%% Dependencies

% add these two code directories to your path
restoredefaultpath;
root_directory = '/Users/svnh2/Desktop/projects';
addpath([root_directory '/general-audio-code']);
addpath(genpath([root_directory '/general-analysis-code']));
addpath([root_directory '/export_fig_v3']);

%%

% Create the TCI sequences.
% Note order is random and so not
% the same as the manuscript.
% S contains info about the
% segments, organized as a cell
% array with one element per segment
% duration tested.
S = create_stimuli;

%% T


T.stim_labels = {};
T.segs = [];
T.orders = [];
T.segorder = struct();
T.sources = unique(S{1}.source);
for i = 1:length(S) % loop through segment durations
    for j = 1:length(S{i}.fname) % loop through sequences of that duration
        [~,stim_label,~] = fileparts(S{i}.fname{j});
        T.stim_labels = cat(1, T.stim_labels, stim_label);
        T.segs = cat(2, T.segs, S{i}.dur);
        T.orders = cat(2, T.orders, j);
    end
end
T.pools = ones(size(T.segs));

%%
T.scramstim_dur = source_duration * n_sources;
T.sourcestim_dur = source_duration;
T.rampwin = 'hann';
T.rampdur = rampdur_sec;

%% Load data

% D: data matrix, time x stimulus x electrode x repetition
% t: time stamps in seconds
% S: structure with stimulus/segment information
TCI_directory = fileparts(fileparts(which('example_data.m')));
data_file = [TCI_directory '/data/example-electrodes.mat'];
X = load(data_file, 'D', 't', 'S');
chnames = {'STG', 'HG'};

%% Cross-context correlationx

% run analysis
% L.diff_context contains the cross-context correlation
% (i.e. the correlation across different contexts)
% format: time x segment duration x electrode
% L.same_context is the noise ceiling
% (i.e. the correlation across identical contexts)
% format: time x segment duration x electrode
directory_to_save_results = [TCI_directory '/results'];
if ~exist(directory_to_save_results, 'dir'); mkdir(directory_to_save_results); end
L = cross_context_corr(D, t, S, 'chnames', chnames, ...
    'output_directory', directory_to_save_results, ...
    'lag_win', [0, 1], 'plot_figure', true);

%% Model-fitting
    
% perform model-fitting
% integration width of best-fitting window is given by M.best_intper_sec
% integration center is given by: M.best_delay_sec_median

% faster version without boundary model parameter
% with single shape parameter, and only 20 tested integration widths
M = modelfit_cross_context_corr(L, 'plot_figure', true, ...
    'shape', 1, 'intper_range', [1/64, 0.5], 'delay_range', [0, 0.5]);

% slower version with boundary model, multiple shape parameters, and 50
% tested integration widths
% M = modelfit_cross_context_corr(L, 'overwrite', false, 'plot_figure', true, ...
%     'shape', [1,2,3,4,5], 'intper_range', [1/64, 0.5], 'delay_range', [0, 0.5], ...
%     'lossfn', 'sqerr', 'nintper', 50, 'boundstrength', [0, 0.25, 0.5, 1, 2]);