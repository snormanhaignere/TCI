function wr = weighted_pearson_corr(X,Y,W)

% Calculates weighted Pearson correlation between corresponding columns of
% X and Y, weighted by W. X, Y, and W should all be the same size.
% 
% 2019-11-18: Created, Sam NH
%
% -- Example -- 
% N = 1000;
% s = randn(N,1);
% w = rand(N,1)*5;
% x = s + randn(N,1).*w;
% y = s + randn(N,1).*w;
% corr(x,y)
% weighted_pearson_corr(x,y,ones(size(x)))
% weighted_pearson_corr(x,y,1./w)

assert(all(size(X)==size(Y)));
assert(all(size(X)==size(W)));

% remove NaNs
xi = isnan(X) | isnan(Y) | isnan(W);
X(xi) = NaN;
Y(xi) = NaN;
W(xi) = 0;

% normalize weights
W = bsxfun(@times, W, 1./sum(W,1));

% weighted means and variances
mx = nansum(W.*X,1);
my = nansum(W.*Y,1);
vx = nansum(W.*(X-mx).^2,1);
vy = nansum(W.*(Y-my).^2,1);

% % weighted covariance
cxy = nansum(W.*(X-mx).*(Y-my));
    
% Pearson correlation
wr = cxy ./ sqrt(vx.*vy);