function F_shifted = add_delays(F, delays, varargin)

% Creates an array with delayed copies of the columns of an input array. Useful
% for convolution and de-convolution.
% 
% 2017-06-06: Created, Sam NH
% 
% -- Very simple example -- 
% 
% imagesc(add_delays(randn(10,1), 0:5));

I.pad_value = 0;
I.direction = 'forward';
I.preserve_singleton = false;
I = parse_optInputs_keyvalue(varargin, I);

% unwrap all dimensions of the feature array but time (the first dimension)
F_dims = size(F);
F = reshape(F, F_dims(1), prod(F_dims(2:end)));

% add delays
% time x features x delays
n_delays = length(delays);
n_tps = F_dims(1);
F_shifted = I.pad_value * ones([F_dims(1), prod(F_dims(2:end)), n_delays]);
original_tps = 1:n_tps;
for i = 1:n_delays
    
    % shifted time points
    switch I.direction
        case 'forward'
            shifted_tps = original_tps + delays(i);
        case 'backward'
            shifted_tps = original_tps - delays(i);
        otherwise
            error('Direction must be forward or backward not %s', I.direction);
    end
    
    % select valid time points
    xi = original_tps < 1 | original_tps > n_tps | shifted_tps < 1 | shifted_tps > n_tps;
    F_shifted( shifted_tps(~xi), :, i ) = F( original_tps(~xi), : );

end

% reshape to original dimensions plus delay
if F_dims(end)==1 && ~I.preserve_singleton
    F_dims = F_dims(1:end-1);
end
F_shifted = reshape(F_shifted, [F_dims, n_delays]);
