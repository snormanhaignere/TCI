function [X1_swapped, X2_swapped] = interleave_oddeven(X1, X2)

% Interleaves odd and even rows of matrices X1 and X2
% 
% 2019-12-05: Created, Sam NH
% 
% -- Example --
% X1 = reshape(1:6, 3, 2);
% X2 = reshape(7:12, 3, 2);
% [X1_swapped, X2_swapped] = interleave_oddeven(X1,X2)

dims = size(X1);

% reshape to matrix
new_dims = [dims(1), prod(dims(2:end))];
X1 = reshape(X1, new_dims);
X2 = reshape(X2, new_dims);

% swap
odd_rows = 1:2:dims(1);
even_rows = 2:2:dims(1);
X1_swapped = nan(new_dims);
X2_swapped = nan(new_dims);
X1_swapped(odd_rows,:) = X1(odd_rows,:);
X1_swapped(even_rows,:) = X2(even_rows,:);
X2_swapped(odd_rows,:) = X2(odd_rows,:);
X2_swapped(even_rows,:) = X1(even_rows,:);

% reshape back
X1_swapped = reshape(X1_swapped, dims);
X2_swapped = reshape(X2_swapped, dims);