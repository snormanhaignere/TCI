% Code that checks model window is correctly parametrized

clc;
intper = 0.2
delay = 0.3
shape = 2;

sr = 100000;
intervalmass = 0.75;
delaypoint = 'median';
tsec = 0:1/sr:2;
% tsec = [];
[h,t_sec,causal] = modelwin('gamma', intper, delay, 'sr', sr, 'delaypoint', delaypoint, ...
    'shape', shape, 'intervalmass', intervalmass, 'tsec', tsec);
plot(t_sec, h);

% measure the integration period
measured_cdf = cumsum(h);
lowtail = 0:0.001:(1-intervalmass)-0.01;
interval = nan(1,length(lowtail));
for i = 1:length(lowtail)
    [~,xi] = unique(measured_cdf);
    interval(i) = diff(interp1(measured_cdf(xi), t_sec(xi), [lowtail(i), lowtail(i) + intervalmass]));
end
clear xi;
measured_intper = min(interval)

% measure the delay
switch delaypoint
    case 'peak'
        [~,xi] = max(h);
        measured_delay = t_sec(xi)
    case 'median'
        measured_cdf = cumsum(h);
        [~,xi] = unique(measured_cdf);
        measured_delay = interp1(measured_cdf(xi), t_sec(xi), 0.5)
        clear xi;
    otherwise
        error('No matching delay point');
end