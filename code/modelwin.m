function [h,t_sec,causal,min_peak] = modelwin(distr, intper_sec, delay_sec, varargin)

% Parametric window with a given integration period and delay.
% 
% Inputs:
% 
% distr: 'gamma' or 'gauss' (determines the shape of the distribution, for
% causal windows, gamma is more appropriate than gauss)
% 
% intper_sec: integration period of the window in seconds (width of the
% central 95% interval)
% 
% delay_sec: delay of the window in seconds, by default the delay
% corresponds to the peak of the window, but this can be changed (see
% optional inputs at the top of the code)
% 
% Outputs: 
% 
% h: the amplitude values for the window
% 
% t_sec: corresponding timestamps
% 
% causal: indicates whether the window is causal or not (i.e. has non-zeros
% values before the 0 timepoint)
% 
% -- Example -- 
% 
% intper_sec = 0.1;
% delay_sec = 0.1;
% distr = 'gamma';
% shape_param_for_gamma = 3; % [1, inf], lower -> more exponential, higher -> more gaussian
% [h,t_sec,causal] = modelwin(distr, intper_sec, delay_sec, ...
%     'plot', true, 'shape', shape_param_for_gamma, 'sr', 1000);
% 
% 2019-11-12: Commented, Sam NH

%% Optional inputs

% default window and sampling rate
I.win = [-intper_sec*3, intper_sec*3] + delay_sec;
I.sr = 100;

% alternatively you can just specify a vector of timepoints
I.tsec = [];

% shape of a Gamma distribution if used
I.shape = 1;

% can optionally zero timepoints before sound onset
I.forcecausal = false;  

% choose how the delay of the window is calculated (only relevant for
% non-Gaussian windows), options are:
% peak: delay is the peak
% median: delay is the median value
% start: delay is the first non-zero point
I.delaypoint = 'peak';

% density interval used to calculate integration period
% intervaltype:
% 'highest': highest density interval (default, the minimal interval of a given mass)
% 'center': central interval, i.e. from [0+(1-intervalmass)/2, 1-(1-intervalmass)/2]
% 'start': starting interval, i.e. from [0 to intervalmass]
I.intervaltype = 'highest';
I.intervalmass = 0.75;

% can optionally return the CDF of the window
I.cdf = false;

% can optionally plot the window
I.plot = false;

% overwrite default parameters with those specified
[I,C] = parse_optInputs_keyvalue(varargin, I);

% set time vector or window depending upon input
if C.tsec % if time-vector is specified, adjust window/sampling rate
    t_sec = I.tsec;
    I.win = I.tsec([1 end]);
    I.sr = 1/diff(I.tsec(1:2));
    assert(all(abs(diff(I.tsec)-diff(I.tsec(1:2)))<1e-6));
else % if time-vector is not specified, determine based on window/sampling rate
    t_sec = (round(I.win(1)*I.sr):round(I.win(2)*I.sr))/I.sr;
end

%% Window

switch distr
    case 'gauss'
        
        error('Mass interval needs to be implemented');
        if ~(strcmp(I.delaypoint, 'peak'))
            error('For Gaussian delaypoint must be "peak"');
        end
        
        % gaussian window
        default_intper = norminv(mass_interval(2),0,1) - norminv(mass_interval(1),0,1);
        sig_sec = intper_sec / default_intper;
        if I.cdf
            h = normcdf(t_sec - delay_sec, 0, sig_sec);
        else
            h = normpdf(t_sec - delay_sec, 0, sig_sec);
        end
        if I.forcecausal
            assert(~I.cdf);
            h(t_sec < 0) = 0;
        end
        h = h / sum(h);
        
        % Gaussian window is never causal
        causal = false;
        min_peak = [];
        
    case 'gamma'
        
        [h, causal] = gamma_reparam(intper_sec, delay_sec, I.shape, ...
            'delaypoint', I.delaypoint, 'cdf', I.cdf, 'intervaltype', I.intervaltype, ...
            'intervalmass', I.intervalmass, 'tsec', t_sec);
        
        if I.forcecausal
            assert(~I.cdf)
            h(t_sec < 0) = 0;
        end
        if ~I.cdf
            h = h / sum(h);
        end
        
    otherwise
        error('No matching distribution');
end

%% Plot

if I.plot
    figure;
    plot(t_sec, h);
    xlim(I.win);
end