function [h_relpower, t_sec, h, causal] = win_power_ratio(...
    segdur_sec, distr, intper_sec, delay_sec, varargin)

% Calculates the expected correlation for TCI correlation analyses by
% measuring the variance a given segment should contribute relative to the
% other segments.
% 
% 2019-13-11: Created, Sam NH

% window and sampling rate
% note the window specifies the window of the output
% which can / should be different from the window used for convolution
I.win = [];
I.sr = NaN;

% alternatively you can just specify a vector of timepoints
% again this is the vector of timepoints of the output
% which can / should be different from the time-points used for convolution
I.tsec = [];

% this is the target sampling rate used for internal computations
I.target_sr = 1000;

% parameters of the window see modelwin.m
I.shape = 1;
I.delaypoint = 'peak';
I.forcecausal = false;
I.centralinterval = 0.75;

% window applied to beginning / end of segment
% see winoverlap.m
I.rampwin = 'none';
I.rampdur = 0.03125; % duration of 

% can optionally plot result
I.plot = false;

% overwrite default parameters with those specified
[I,C] = parse_optInputs_keyvalue(varargin, I);

%%  compute window

% ensure there are an integer number of samples in the segment
internal_sr = ceil(segdur_sec * I.target_sr) / segdur_sec;
segdur_smps = checkint(segdur_sec * internal_sr);

% create the window
[h, t, causal] = winoverlap(segdur_sec, distr, intper_sec, delay_sec, ...
     'shape', I.shape, 'forcecausal', I.forcecausal, 'delaypoint', I.delaypoint, ...
     'plot', false,'internal_sr', internal_sr, 'rampwin', I.rampwin, 'rampdur', I.rampdur, ...
     'centralinterval', I.centralinterval);
h(h<0) = 0;

% create delayed copies
N_delays = ceil(length(h)/segdur_smps);
shifts = (-N_delays:N_delays)*segdur_smps;
h_delayed = add_delays(h, shifts);

% normalize by summed power across filters
h_relpower = bsxfun(@times, h.^2, 1./sum(h_delayed.^2,2));

%% Interpolate

if ~C.sr && ~C.win && ~C.tsec
    t_sec = t;
else
    % set time vectory or window depending upon input
    if C.tsec % if time-vector is specified, adjust window/sampling rate
        t_sec = I.tsec;
        I.win = I.tsec([1 end]);
        I.sr = 1/diff(I.tsec(1:2));
        assert(all(abs(diff(I.tsec)-diff(I.tsec(1:2)))<1e-6));
    else % if time-vector is not specified, determine based on window/sampling rate
        t_sec = (round(I.win(1)*I.sr):round(I.win(2)*I.sr))/I.sr;
    end
    if any(isnan(h_relpower(:))) || any(isnan(h(:))) || any(isnan(h_delayed(:)))
        error('NaN values');
    end
    h_relpower = interp1(t, h_relpower, t_sec(:), 'pchip');
    h = interp1(t, h, t_sec(:), 'pchip');
    xi = t_sec<t(1) | t_sec>t(end);
    h_relpower(xi) = 0;
    h(xi) = 0;
end

%% Plot

if I.plot
    figure;
    subplot(2,1,1);
    h = plot(t_sec, [h, h_relpower], 'LineWidth', 2);
    xlim(t_sec([1,end]));
    legend(h, {'Orig', 'Rel Power'}, 'Location', 'Best');
    subplot(2,1,2);
    plot(t, h_delayed, 'LineWidth', 2);
    hold on;
    plot(t, sum(h_delayed,2), 'k--', 'LineWidth', 2);
    xlim(t([1,end]));
    title('All win');
end