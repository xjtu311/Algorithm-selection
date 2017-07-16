function logrank_stat = logrank_statistic(N,O,N_1,O_1)
%=== Compute logrank statistic for the specified "numbers at risk" (N),
%observed events (O), "numbers at risk" in group 1 (N_1), and observed
%events in group 1 (O_1)
%
% If N(i) <= 1, we don't include it in the sum (undefined variance, 0/0)

E = O .* (N_1+0.0) ./ N;

zero_idx = find(N <= 1);
nonzero_idx = find(N > 1);
V(nonzero_idx) = E(nonzero_idx) .* (1 - ((N_1(nonzero_idx)+0.0) ./ N(nonzero_idx))) .* (N(nonzero_idx)-O(nonzero_idx)) ./ (N(nonzero_idx)-1.0);
V(zero_idx) = 0;

denominator = sqrt(sum(V(nonzero_idx)));
numerator = sum(O_1(nonzero_idx)-E(nonzero_idx));

if denominator == 0
    if numerator == 0
        logrank_stat = 0;
    else
        error 'Division by zero in function logrank_statistic'
    end
else
    logrank_stat = numerator/denominator;
end
assert(~isnan(logrank_stat));