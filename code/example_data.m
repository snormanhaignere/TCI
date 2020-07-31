% Demonstrates code
% 
% 2019-11-26: Created, Sam NH

%% Load data

TCI_directory = fileparts(fileparts(which('example_data.m')));
data_file = [TCI_directory '/data/example-electrodes.mat'];
load(data_file, 'D', 't', 'S', 'chnames');

%% Cross-context correlation

clc;
tic;
root_directory = '/Users/svnh2/Desktop/projects';
addpath([root_directory '/general-analysis-code']);
addpath([root_directory '/export_fig_v3']);
directory_to_save_results = [TCI_directory '/results/v6'];
L = cross_context_corr(D, t, S, 'chnames', chnames, ...
    'output_directory', directory_to_save_results, ...
    'interleave_diffdur', false, 'boundary', 'any', ...
    'overwrite', true, 'plot_figure', true, 'lag_win', [0, 1], ...
    'excludewin', []);
toc;

%% Model-fitting
    
clc;
M = modelfit_cross_context_corr(L, 'overwrite', false, 'plot_figure', true, 'plot_jack', true, ...
    'shape', 1, 'lossfn', 'unbiased-sqerr', 'nullsmps', 10, 'intper_range', [1/16, 1], ...
    'skipnullpreds', true, 'delay_range', [0, 0.1], 'skipsmppreds',  false);
