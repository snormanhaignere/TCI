function shuff = pseudorand(N)

if N == 1
    shuff = 1;
    return;
end

while true
    shuff = Shuffle(1:N);
    prev1 = 0:N-1;
    prev2 = nan(1, N);
    prev2(shuff(1:end)) = [0, shuff(1:end-1)];
    if all(prev1 ~= prev2)
        break;
    end
end