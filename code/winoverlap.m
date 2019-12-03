function [overlap,t_sec,causal] = winoverlap(segdur_sec, distr, intper_sec, delay_sec, varargin)

% Convolves a model window (see modelwin.m) with a segment of a given
% duration. Useful for modeling TCI data.
% 
% Inputs:
% 
% segdur_sec: duration of the segment
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
% overlap: result of convolution with segment
% 
% t_sec: corresponding timestamps
% 
% causal: returns whether the window was causal
% 
% -- Example --
% 
% segdur_sec = 0.2;
% distr = 'gamma';
% intper_sec = 0.1;
% delay_sec = 0.1;
% shape_param_for_gamma = 3; % [1, inf], lower -> more exponential, higher -> more gaussian
% winoverlap(segdur_sec, distr, intper_sec, delay_sec, ...
%     'plot', true, 'shape', shape_param_for_gamma);
% 
% 2019-11-13: Cleaned up, created internal time vector, Sam NH

%% Optional inputs

time_boundary = [-1,1]*(3*segdur_sec + 4*intper_sec + abs(delay_sec));
% modelwin_time_boundary = [-intper_sec*3, intper_sec*3]+delay_sec;
% time_boundary = [min(seg_time_boundary(1),modelwin_time_boundary(1)), max(seg_time_boundary(2), modelwin_time_boundary(2))];

% window and sampling rate
% note the window specifies the window of the output
% which can / should be different from the window used for convolution
I.win = time_boundary;
I.sr = 1000;

% alternatively you can just specify a vector of timepoints
% again this is the vector of timepoints of the output
% which can / should be different from the time-points used for convolution
I.tsec = [];

% internal sampling rate used for convolution
I.internal_sr = 1000;

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

% window applied to beginning / end of segment
I.rampwin = 'none';
I.rampdur = 0.03125; % duration of 

% can optionally plot the result
I.plot = false;

% overwrite default parameters with those specified
[I,C] = parse_optInputs_keyvalue(varargin, I);

%% Perform convolution

% internal time vectory used for convolution
internal_tsec = time_boundary(1):1/I.internal_sr:time_boundary(2);

% adjust vector so there is a zero timepoint
[~,zero_tp] = min(abs(internal_tsec-0));
internal_tsec = internal_tsec - internal_tsec(zero_tp);

% rectangular segment
switch I.rampwin
    case 'none'
        seg = zeros(size(internal_tsec));
        xi = internal_tsec <= segdur_sec & internal_tsec >= 0;
        seg(xi) = 1;

    case 'hann'
        seg = zeros(size(internal_tsec));
        xi1 = find(internal_tsec <= segdur_sec & internal_tsec >= 0, 1, 'first');
        xi2 = find(internal_tsec <= segdur_sec & internal_tsec >= 0, 1, 'last');
        rampdur_smp = round(I.rampdur*I.sr);
        halframpdur_smp1 = floor(rampdur_smp/2);
        halframpdur_smp2 = rampdur_smp - halframpdur_smp1;
        xi = xi1-halframpdur_smp1:xi2+halframpdur_smp2;
        seg(xi) = ramp_hann(ones(length(xi),1), I.rampdur, I.sr);
        
        %     case 'hann'
        %         seg = zeros(size(internal_tsec));
        %         xi = internal_tsec <= segdur_sec & internal_tsec >= 0;
        %         seg(xi) = ramp_hann(ones(sum(xi),1), I.rampdur, I.sr);
        %
        %     case 'hannextend'
        %         seg = zeros(size(internal_tsec));
        %         xi = internal_tsec <= segdur_sec & internal_tsec >= 0;
        %         halframpdur_smp = round(I.rampdur*I.sr/2);
        %         x = ramp_hann(ones(sum(xi)+halframpdur_smp*2,1), I.rampdur, I.sr);
        %         seg(xi) = x(1+halframpdur_smp:end-halframpdur_smp);


    otherwise
        error('No matching value for rampwin');
end

% model window
[h,~,causal] = modelwin(distr, intper_sec, delay_sec,...
    'shape', I.shape, 'forcecausal', I.forcecausal, 'delaypoint', I.delaypoint, ...
    'tsec', internal_tsec, 'intervalmass', I.intervalmass, ...
    'intervaltype', I.intervaltype);

% convolve
internal_overlap = myconv(seg', h', 'causal', false, 'central_point', zero_tp);

%% Interpolate

if ~C.sr && ~C.win && ~C.tsec
    overlap = internal_overlap;
    t_sec = internal_tsec;
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
    overlap = interp1(internal_tsec, internal_overlap, t_sec, 'pchip');
end

%% Plot

if I.plot
    figure;
    plot(t_sec, overlap);
    xlim([t_sec(1), t_sec(end)])
end