clc

n_reps = 1000;

r = nan(n_reps,1);
r1 = nan(n_reps,1);
r2 = nan(n_reps,1);
for i = 1:n_reps
    N = 100;
    noisestrength = 1;
    same_context = rand(N,1);
    same_context_noisy = same_context + randn(N,1)*noisestrength;
    true_model = rand(N,1);
    diff_context_noisy = same_context.*true_model + randn(N,1)*noisestrength;
    
    p = true_model.*same_context_noisy;
    p1 = 0.8*true_model.*same_context_noisy;
    p2 = 1.2*true_model.*same_context_noisy;
   
    r(i) = 1-normalized_squared_error(p,diff_context_noisy);
    r1(i) = 1-normalized_squared_error(p1,diff_context_noisy);
    r2(i) = 1-normalized_squared_error(p2,diff_context_noisy);
    r(i) = sum((p-diff_context_noisy).^2);
    r1(i) = sum((p1-diff_context_noisy).^2);
    r2(i) = sum((p2-diff_context_noisy).^2);
end

mean(r)
mean(r1)
mean(r2)