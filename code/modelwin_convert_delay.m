function converted_delay_sec = modelwin_convert_delay(...
    intper_sec, delay_sec, shape, delaypoint, varargin)

% Calculates the peak or median delay of a Gamma distributed window with a
% given integration period, shape and "delay" where delay is determined by
% the first non-zero point. This function is useful because the
% modelfitting functions asumme the start is at the first non-zero point.
% 
% -- Example --
% 
% intper_sec = 1;
% delay_sec_start = 0;
% shape = 1;
% delaypoint = 'median';
% converted_delay_sec = modelwin_convert_delay(intper_sec, delay_sec_start, shape, delaypoint);
% 
% % plot delay point
% [h, t] = modelwin('gamma', intper_sec, delay_sec_start, 'shape', shape, 'delaypoint', 'start');
% plot(t, h, 'LineWidth', 2);
% hold on;
% yL = ylim;
% lineh = plot([1,1]*converted_delay_sec, yL, 'r--', 'LineWidth', 2);
% legend(lineh, delaypoint);
% 
% 2019-11-26: Created, Sam NH

clear I;

% density interval used to calculate integration period
% intervaltype:
% 'highest': highest density interval (default, the minimal interval of a given mass)
% 'center': central interval, i.e. from [0+(1-intervalmass)/2, 1-(1-intervalmass)/2]
% 'start': starting interval, i.e. from [0 to intervalmass]
I.intervaltype = 'highest';
I.intervalmass = 0.75;
I = parse_optInputs_keyvalue(varargin, I);

% set interval
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
        mass_interval = [low_tail, high_tail];
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
switch delaypoint
    case 'peak'
        converted_delay_sec = max((a-1)*b,0)*r + delay_sec;
    case 'median'
        converted_delay_sec = gaminv(0.5,a,b)*r + delay_sec;
    otherwise
        error('delaypoint %s is not valid\n', I.delaypoint);
end