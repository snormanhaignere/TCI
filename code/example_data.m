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
directory_to_save_results = [TCI_directory '/results/v4'];
L = cross_context_corr_further_simplified(D, t, S, 'chnames', chnames, ...
    'output_directory', directory_to_save_results, 'boundary', 'any', ...
    'overwrite', false, 'plot_figure', true, 'interleave', true, 'lag_win', [0, 1]);
toc;

%% Model-fitting

M = modelfit_cross_context_corr(L, 'overwrite', true, ...
    'shape', [1,2,3,5], 'lossfn', 'nse');