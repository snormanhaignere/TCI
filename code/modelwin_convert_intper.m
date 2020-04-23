function converted_intper_sec = modelwin_convert_intper(...
    intper_sec, shape, orig_intervalmass, new_intervalmass, varargin)

% Given a window calculates the equivalent integration period 
% for a new mass.
% 
% 2019-03-12: Created, Sam NH

clear I;

% density interval used to calculate integration period
% intervaltype:
% 'highest': highest density interval (default, the minimal interval of a given mass)
% 'center': central interval, i.e. from [0+(1-intervalmass)/2, 1-(1-intervalmass)/2]
% 'start': starting interval, i.e. from [0 to intervalmass]
I.intervaltype = 'highest';
I = parse_optInputs_keyvalue(varargin, I);

% set interval
switch I.intervaltype
    case 'center'
        low_tail = (1-orig_intervalmass)/2;
        high_tail = 1-low_tail;
        orig_mass_interval = [low_tail, high_tail];
        low_tail = (1-new_intervalmass)/2;
        high_tail = 1-low_tail;
        new_mass_interval = [low_tail, high_tail];
    case 'start'
        low_tail = 0;
        high_tail = orig_intervalmass;
        orig_mass_interval = [low_tail, high_tail];
        low_tail = 0;
        high_tail = new_intervalmass;
        new_mass_interval = [low_tail, high_tail];
    case 'highest'
        load('highest_density_interval_gamma.mat', 'H');
        xi = abs(H.mass-orig_intervalmass)<1e-6;
        yi = abs(H.shape-shape)<1e-6;
        assert(sum(xi)==1);
        assert(sum(yi)==1);
        low_tail = H.lowtail(xi,yi);
        high_tail = low_tail + orig_intervalmass;
        orig_mass_interval = [low_tail, high_tail]; % Original CDF interval
        orig_physical_interval = H.interval(xi,yi);
        xi = abs(H.mass-new_intervalmass)<1e-6;
        assert(sum(xi)==1);
        low_tail = H.lowtail(xi,yi);
        high_tail = low_tail + new_intervalmass;
        new_mass_interval = [low_tail, high_tail]; % New CDF Interval
        new_physical_interval = H.interval(xi,yi);
    otherwise
        error('No matching interval type');
end

% gaussian window
% sigma corresponding to 95%
a = shape;
b = 1/shape;

% ratio which to scale stimulus
orig_default_intper = gaminv(orig_mass_interval(2),a,b) - gaminv(orig_mass_interval(1),a,b);
assert(abs(orig_physical_interval - orig_default_intper)<1e-6);
new_default_intper = gaminv(new_mass_interval(2),a,b) - gaminv(new_mass_interval(1),a,b);
assert(abs(new_physical_interval - new_default_intper)<1e-6);
converted_intper_sec = intper_sec * new_default_intper/orig_default_intper;