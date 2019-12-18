function M_new = split_dim_into_fields(M, f_orig, dim, split1, split2, f_new1, f_new2, split1_rs, split2_rs, squeeze_split1, squeeze_split2)

% Helper function for modelfit_cross_context_corr.m
% 
% 2019-12-17: Created Sam NH



% checks
assert(length(f_orig)==length(f_new1));
assert(length(f_orig)==length(f_new2));

M_new = M;
for i = 1:length(f_orig)
    d = size(M.(f_orig{i}));
    if nargin < 8 || isempty(split1_rs)
        d_rs1 = [d(1:dim-1), length(split1), d(dim+1:end)]; 
    else
        d_rs1 = [d(1:dim-1), split1_rs, d(dim+1:end)];
    end
    if nargin < 9 || isempty(split1_rs)
        d_rs2 = [d(1:dim-1), length(split2), d(dim+1:end)]; 
    else
        d_rs2 = [d(1:dim-1), split2_rs, d(dim+1:end)];
    end
    M_new.(f_new1{i}) = reshape(index(M.(f_orig{i}), dim, split1), d_rs1);
    M_new.(f_new2{i}) = reshape(index(M.(f_orig{i}), dim, split2), d_rs2);
    if nargin >= 10 && squeeze_split1 && length(split1)==1
        M_new.(f_new1{i}) = squeeze_dims(M_new.(f_new1{i}), dim);
    end
    if nargin >= 11 && squeeze_split2 && length(split2)==1
        M_new.(f_new2{i}) = squeeze_dims(M_new.(f_new2{i}), dim);
    end
end
clear M;