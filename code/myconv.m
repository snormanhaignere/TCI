function y = myconv(sig, h, varargin)

% Similar to conv.m but has default behaviors that I find more intuitive. By
% default the function assumes that the kernel (h) is causal, and thus is
% flipped slid against the signal vector (sig), at each step computing the
% dot-product. The signal is zero-padded with sufficient padding to avoid
% wrap-around effects. The padded segments are then removed after convolution.
% The signal can be a matrix, in which case the 1D kernel is convolved
% separately with each column of the signal.
% 
% To make the kernel non-causal, add the optional arguments:
% 
% myconv(..., 'causal', true).
% 
% In this case the kernel is shifted so that the "central timepoint" is the
% first entry of the kernel, and the kernel is then convolved with signal as per
% usual (flipped and slid against the signal). By default the central timpoint
% is the center point of an odd signal (i.e. 3rd element of a signal of length
% 5), or the element just to the right of center for an even signal
% (length/2+1). The central timepoint can be manually specified by including the
% optional arguments:
% 
% myconv(...,'central_timepoint', DESIREDTIMEPOINT).
% 
% It is also possible to normalize the output by the sum of the kernel being
% multiplied by the non-padded signal. To do this, add the optional arguments:
% 
% myconv(...,'norm', true)
% 
% myconv(sig, h, 'causal', false) is equivalent to conv(sig, h, 'same')
% 
% 2016-09-25: Created, Sam NH
% 
% 2020-06-17: Fixed a major bug (fft was fft2!)

%% Optional arguments

% whether or not the kernel h is assumed to be causal
I.causal = true;

% whether to normalize by the sum of kernel multiplying the nonpadded signal
I.norm = false;

% central point for non-causal filters
if mod(length(h),2)==1 % odd length kernel -> central point
    I.central_point = ceil(length(h)/2);
else % even length kernel -> just right of center
    I.central_point = length(h)/2 + 1;
end

% parse changes to the above optional arguments
I = parse_optInputs_keyvalue(varargin, I);
    
%% Padding

assert(isvector(h));

% pad signal with zeros
sig_pad = [zeros(size(h,1), size(sig,2)); sig];

% pad kernel to match size of padded signal
h_pad = [h; zeros(size(sig_pad,1)-size(h,1),1)];

% shift padded kernel so the central timepoint is the first element
if ~I.causal
    h_pad = circshift(h_pad, -(I.central_point-1));
end

%% Convolution and normalization
    
% convolve in the frequency domain
y = ifft(fft(sig_pad) .* (fft(h_pad) * ones(1,size(sig_pad,2))));

% remove padded segments
y = y(size(h,1)+1:end,:);

% normalize by the sum of kernel multiplying the nonpadded signal
if I.norm
    xi = optInputs(varargin, 'norm');
    assert(length(xi)==1);
    varargin = varargin(setdiff(1:length(varargin), [xi,xi+1]));
    y_ones = myconv(ones(size(sig)), h, 'norm', false, varargin{:});
    y = y./y_ones;
end
    