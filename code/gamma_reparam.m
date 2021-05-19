function [h, causal, a, b, r, c, t_sec, min_peak] = gamma_reparam(intper_sec, delay_sec, shape, varargin)

% Reparametrize gamma distribution in terms of 
% integration period, delay, and shape

I.delaypoint = 'peak';
I.intervaltype = 'highest';
I.intervalmass = 0.75;
I.cdf = false;
I.win = [-intper_sec*3, intper_sec*3] + delay_sec;
I.sr = 100;
I.tsec = [];
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

% CDF interval corresponding to desired cumulative mass
% i.e. 10-60% or 20-70%
switch I.intervaltype
    case 'center'
        low_tail = (1-I.intervalmass)/2;
        high_tail = 1-low_tail;
        mass_interval = [low_tail, high_tail];
    case 'start'
        low_tail = 0;
        high_tail = I.intervalmass;
        mass_interval = [low_tail, high_tail];
    case 'highest'
        load('highest_density_interval_gamma.mat', 'H');
        xi = abs(H.mass-I.intervalmass)<1e-6;
        yi = abs(H.shape-shape)<1e-6;
        assert(sum(xi)==1);
        assert(sum(yi)==1);
        low_tail = H.lowtail(xi,yi);
        high_tail = low_tail + I.intervalmass;
        mass_interval = [low_tail, high_tail]; % CDF interval
    otherwise
        error('No matching interval type');
end

% gaussian window
% sigma corresponding to 95%
a = shape;
b = 1/shape;

% ratio which to scale stimulus
default_intper = gaminv(mass_interval(2),a,b) - gaminv(mass_interval(1),a,b);
r = intper_sec/default_intper;

% offset to adust delay
switch I.delaypoint
    case 'peak'
        min_peak = max((a-1)*b,0)*r;
    case 'median'
        min_peak = gaminv(0.5,a,b)*r;
    case 'start'
        min_peak = 0;
    otherwise
        error('delaypoint %s is not valid\n', I.delaypoint);
end
c = delay_sec - min_peak;
if delay_sec < min_peak
    causal = false;
else
    causal = true;
end

% gamma distribution
if I.cdf
    h = gamcdf((t_sec-c)/r, a, b);
else
    h = gampdf((t_sec-c)/r, a, b);
end

