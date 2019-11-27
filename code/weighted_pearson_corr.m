function wr = weighted_pearson_corr(X,Y,W)

% Calculates weighted Pearson correlation between corresponding columns of
% X and Y, weighted by W.
% 
% 2019-11-18: Created, Sam NH
% 
% 2019-11-27: Got rid of some house-keeping that was slowing thiings down.
% Also using bsxfun for matrix expansion.
%
% -- Example -- 
% N = 1000;
% s = randn(N,1);
% w = rand(N,1)*5;
% x = s + randn(N,1).*w;
% y = s + randn(N,1).*w;
% invW = 1./w;
% corr(x,y)
% weighted_pearson_corr(x,y,ones(size(x))/N)
% weighted_pearson_corr(x,y,invW/sum(invW))

% weighted means
Mx = sum(bsxfun(@times, W, X),1);
My = sum(bsxfun(@times, W, Y),1);

% demeaned variables
X_demean = bsxfun(@minus, X, Mx);
Y_demean = bsxfun(@minus, Y, My);

% weighted variances
vx = sum(bsxfun(@times, W, X_demean.^2),1);
vy = sum(bsxfun(@times, W, Y_demean.^2),1);

% % weighted covariances
XY = bsxfun(@times, X_demean, Y_demean);
cxy = sum(bsxfun(@times, W, XY), 1);
    
% Pearson correlation
wr = cxy ./ sqrt(bsxfun(@times, vx, vy));