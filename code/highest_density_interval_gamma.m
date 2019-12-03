% Calculates and saves the highest density interval for a Gamma
% distribution of a given shape and mass
% 
% Used in modelwin.m for defining integraitoin periods.
% 
% 2019-12-3: Created, Sam NH

% the total mass in the interval
H.mass = 0.05:0.05:0.95;

% shape of the Gamma distrubiton
H.shape = 1:0.5:10;

H.lowtail = nan(length(H.mass), length(H.shape));
H.interval = nan(length(H.mass), length(H.shape));
for i = 1:length(H.mass)
    for j = 1:length(H.shape)
        
        shape = H.shape(j);
        mass = H.mass(i);
        a = shape;
        b = 1/shape;
        
        % the lower and upper tail
        low_tail = 0:0.001:mass-0.001;
        high_tail = low_tail + mass;
        
        % the interval corresponding for a given lower/upper tail of fixed mas
        interval = gaminv(high_tail,a,b) - gaminv(low_tail,a,b);
        
        % find and store minimum interval
        [~,xi] = min(interval);
        H.lowtail(i,j) = low_tail(xi);
        H.interval(i,j) = interval(xi);
        
    end
end

% save results
save('highest_density_interval_gamma.mat', 'H');