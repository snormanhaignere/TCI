function M_new = split_dim_into_fields(M, f_orig, dim, splits, f_new, varargin)

clear I;
I.shapes = {};
I.squeeze = false(0);
I = parse_optInputs_keyvalue(varargin, I);

% original dimensionality
d = size(M.(f_orig));

assert(length(splits)==length(f_new));
assert(all(unique(cat(2, splits{:}))==(1:size(M.(f_orig), dim))));
M_new = M;
for j = 1:length(splits)
    if isempty(I.shapes)
        M_new.(f_new{j}) = index(M.(f_orig), dim, splits{j});
    else
        assert(length(I.shapes)==length(splits));
        M_new.(f_new{j}) = reshape(...
            index(M.(f_orig), dim, splits{j}), ...
            [d(1:dim-1), I.shapes{j}, d(dim+1:end)]);
    end
    if ~isempty(I.squeeze) && I.squeeze(j)
        assert(length(I.squeeze)==length(splits));
        M_new.(f_new{j}) = squeeze_dims(M_new.(f_new{j}), dim);
    end
end