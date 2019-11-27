function wnse = weighted_nse(X,Y,W)

% Calculates weighted NSE between corresponding columns of
% X and Y, weighted by W.
% 
% 2019-11-18: Created, Sam NH
% 
% 2019-11-27: Fixed bug (using variance not power), and got rid of some
% house-keeping that was slowing thiings down. Also using bsxfun for matrix
% expansion.
%
% -- Example -- 
% N = 1000;
% s = randn(N,1);
% w = rand(N,1)*5;
% x = s + randn(N,1).*w;
% y = s + randn(N,1).*w;
% invW = 1./w;
% 1-normalized_squared_error(x,y)
% weighted_nse(x,y,ones(size(x))/N)
% weighted_nse(x,y,invW/sum(invW))

% weighted means
Mx = sum(bsxfun(@times, W, X),1);
My = sum(bsxfun(@times, W, Y),1);

% weighted power and cross power
Px = sum(bsxfun(@times, W, X.^2),1);
Py = sum(bsxfun(@times, W, Y.^2),1);
Pxy = sum(bsxfun(@times, W, bsxfun(@times, X, Y)),1);

% calculate NSE
a = bsxfun(@plus, Px, Py) - 2 * Pxy;
b = bsxfun(@plus, Px, Py) - 2 * bsxfun(@times, Mx, My);
wnse = a ./ b;