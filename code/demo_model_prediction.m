segdur_sec = 0.8;
intper_sec = 0.1;
delay_sec = 0;
shape = 3; % 1 2 3 5    

sr = 100;
intervalmass = 0.75;
delaypoint = 'start';
tsec = 0:1/sr:2;
[h, t] = win_power_ratio(segdur_sec, 'gamma', intper_sec, delay_sec, 'sr', sr, 'delaypoint', delaypoint, ...
    'shape', shape, 'intervalmass', intervalmass, 'tsec', tsec);

plot(t, h)