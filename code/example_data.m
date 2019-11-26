% Demonstrates code
% 
% 2019-11-26: Created, Sam NH

%% Load data

TCI_directory = fileparts(fileparts(which('example_data.m')));
data_file = [TCI_directory '/data/example-electrodes.mat'];
load(data_file, 'D', 't', 'S', 'chnames');

%% Cross-context correlation

directory_to_save_results = [TCI_directory '/results'];
L = cross_context_corr(D, t, S, 'chnames', chnames, ...
    'output_directory', directory_to_save_results, 'boundary', 'none');

%% Model-fitting

M = modelfit_cross_context_corr_fast(L, 'overwrite', true);