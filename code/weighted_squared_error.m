function we = weighted_squared_error(X,Y,W)

assert(all(size(X)==size(Y)));
assert(all(size(X)==size(W)));

% remove NaNs
xi = isnan(X) | isnan(Y) | isnan(W);
X(xi) = NaN;
Y(xi) = NaN;
W(xi) = 0;

% normalize weights
W = bsxfun(@times, W, 1./sum(W,1));

% weighted squared error
we = nansum(W.*(X - Y).^2);