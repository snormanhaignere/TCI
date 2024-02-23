function y = checkint(x, tol)

% Checks to make sure that the input is an integer within some degree of
% tolerance. Returns the nearest integer if it is within the tolerance range.

if nargin < 2
    tol = 1e-10;
end
y = round(x);
assert(all(abs(x - y) < tol));


