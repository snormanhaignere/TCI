P = 1000;
r = nan(P,1);
for i = 1:P
    N = 100;
    s = randn(N,1);
    y1 = s + randn(N,1)*2;
    y2 = s + randn(N,1)*2;
    r(i) = corr(y1,y2);
end

stderr_from_samples(r)
r2rerr(mean(r), N)

