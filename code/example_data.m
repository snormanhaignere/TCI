% Illustrates how to apply the code to the example electrode responses from 
% this paper:
% 
% Norman-Haignere SV, Long LK, Devinsky O, Doyle W, Irobunda I, Merricks
% EM, Feldstein NA, Schevon CA, Flinker A, Mesgarani N. Multiscale temporal
% integration organizes hierarchical computation in human auditory cortex.
% Nature Human Behavior.

%% Dependencies

% These code directories need to be added to your path
root_directory = '/Users/svnh2/Desktop/projects';

% see github.com/snormanhaignere/general-audio-code
addpath([root_directory '/general-audio-code']);

% see github.com/snormanhaignere/general-analysis-code
addpath(genpath([root_directory '/general-analysis-code']));

% code needed to save plots
% see mathworks.com/matlabcentral/fileexchange/23629-export_fig
% will also need Ghostscript (https://www.ghostscript.com/)
addpath(genpath([root_directory '/export_fig_v3']));

%% Load data for example electrodes

% D: data matrix, time x stimulus x electrode x repetition
% t: time stamps in seconds
% S: structure with stimulus/segment information
% The format of S is unfortunately somewhat opaque (which we're working on).
% Contact Sam Norman-Haignere (samuel_norman-haignere@urmc.rochester.edu)
% if you need assistance. 
TCI_directory = fileparts(fileparts(which('example_data.m')));
data_file = [TCI_directory '/data/example-electrodes.mat'];
load(data_file, 'D', 't', 'S');
chnames = {'HG', 'STG'};

%% Estimate cross-context correlation

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
    'lag_win', [-2, 4], 'plot_figure', true, 'plot_win', [-0.5, 1.5]);

%% Model-fitting using cross-context correlation
    
% Use model to estimate integration window from cross-context correlation
% Integration width of best-fitting window is given by M.best_intper_sec
% Integration center is given by: M.best_delay_sec_median

% Faster version without boundary model parameter
% with single shape parameter, and only 20 tested integration widths
M = modelfit_cross_context_corr(L, 'plot_figure', true, ...
    'shape', 1, 'intper_range', [1/64, 0.5], 'delay_range', [0, 0.5]);  

% Slower version with boundary model, multiple shape parameters, and 100
% tested integration widths over a slightly wider range (0, 1 second)
M = modelfit_cross_context_corr(L, 'overwrite', false, 'plot_figure', false, ...
    'shape', [1,2,3,4,5], 'intper_range', [1/32, 1], 'delay_range', [0, 0.5], ...
    'lossfn', 'unbiased-sqerr', 'nintper', 100, 'boundstrength', [0, 0.25, 0.5, 1, 2]);
for i = 1:2
    fprintf('%s: width=%.0f ms, center=%.0f ms\n', ...
        chnames{i}, M.best_intper_sec(i)*1000, M.best_delay_sec_median(i)*1000);
end