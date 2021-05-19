%% Load data

% D: data matrix, time x stimulus x electrode x repetition
% t: time stamps in seconds
% S: structure with stimulus information
TCI_directory = fileparts(fileparts(which('example_data.m')));
data_file = [TCI_directory '/data/example-electrodes.mat'];
load(data_file, 'D', 't', 'S', 'chnames');

%% Cross-context correlationx

% L.diff_context contains the cross-context correlation
% L.same_context is the noise ceiling
root_directory = '/Users/svnh2/Desktop/projects';
addpath([root_directory '/general-analysis-code']);
addpath([root_directory '/export_fig_v3']);
directory_to_save_results = [TCI_directory '/results/v1'];
L = cross_context_corr(D, t, S, 'chnames', chnames, ...
    'output_directory', directory_to_save_results, ...
    'boundary', 'any', 'lag_win', [0, 1], 'plot_figure', true);

%% Model-fitting
    
% estimated integration width in seconds: M.best_intper_sec
% estimated integration center/delay in seconds: M.best_delay_sec_median
M = modelfit_cross_context_corr(L, 'overwrite', false, 'plot_figure', true, ...
    'shape', [1,2,3,4,5], 'intper_range', [1/16, 1], 'delay_range', [0, 0.5]);