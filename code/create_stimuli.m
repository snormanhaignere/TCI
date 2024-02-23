function S = create_stimuli

% Creates two different randomly ordered segment sequences
% note the order is random and so not the same as in the manuscript
% 
% Returns a structure S with the information about the segments.
% This structure is needed for the analysis.
% 
% Relies heavily on concat_with_crossfade.m

% sounds from which to create segments
% the 10 sounds here are the same as those from
% the NHB paper
TCI_directory = fileparts(fileparts(which('demo_TCI_paradigm.m')));
source_directory = [TCI_directory '/stimuli/source_sounds'];
source_stimuli = mydir(source_directory);
n_sources = length(source_stimuli);

% directory to save sequences to
output_directory = [TCI_directory '/stimuli/segment-sequences'];

% different segment durations to test
segment_durations = [1/32, 1/16, 1/8, 1/4, 1/2, 1, 2];
n_segdur = length(segment_durations);

% number of different random orderings
n_orders = 2;

% going to use the first 2-seconds of every source file, discarding the
% rest if the source happens to be longer
source_duration = 2; 

% sampling rate
sr = 20e3;

% add zero padding to source
% needed for cross-fading with first/last segment
pad_dur = 1/32;

S = cell(1, n_segdur);
for i = 1:n_segdur
    
    % create structure that species the segments to use
    S{i}.dur = segment_durations(i);
    S{i}.source = {};
    S{i}.onset = [];
    S{i}.level = [];
    S{i}.normwin = [];
    S{i}.directory = source_directory;
    n_segs_per_source = checkint(source_duration / segment_durations(i));
    for j = 1:n_sources
        for k = 1:n_segs_per_source
            % file name of the segment
            S{i}.source = cat(1, S{i}.source, source_stimuli{j});
            
            % onset of the segment relative to the source
            S{i}.onset = cat(1, S{i}.onset, (k-1)*segment_durations(i)+pad_dur);
        end
    end
    % RMS normalization is going to be applied to the source not the segments
    % otherwise segment level would vary with segment duration
    S{i}.level = 0.025;
    S{i}.normwin = [0, 2];
    
    % concatenate segments together in random orders
    rampdur_sec = 1/32;
    n_segs = length(S{i}.source);
    S{i}.order = nan(n_segs, n_orders);
    S{i}.fname = cell(1, n_orders);
    for j = 1:n_orders
        
        % create a random segment order
        S{i}.order(:,j) = randperm(n_segs);
        
        % concatenate the segments using the above order
        stim = concat_with_crossfade(S{i}, S{i}.order(:,j), sr, rampdur_sec, 'RemoveVeryLow', false, 'pad', pad_dur);
        
        % write to file
        S{i}.fname{j} = mkpdir([output_directory '/seg-' num2str(segment_durations(i)) '-order' num2str(j) '.wav']);
        audiowrite(S{i}.fname{j}, stim, sr);
    end
end