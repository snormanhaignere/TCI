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
I.centralinterval = 0.75;
I = parse_optInputs_keyvalue(varargin, I);

% calculate central interval
low_tail = (1-I.centralinterval)/2;
high_tail = 1-low_tail;
central_interval = [low_tail, high_tail];

% gaussian window
% sigma corresponding to 95%
a = shape;
b = 1/shape;

% ratio which to scale stimulus
default_intper = gaminv(central_interval(2),a,b) - gaminv(central_interval(1),a,b);
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