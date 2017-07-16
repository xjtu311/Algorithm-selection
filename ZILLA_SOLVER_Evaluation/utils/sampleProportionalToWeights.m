function idxs = sampleProportionalToWeights(weights, N)
%=== Sample N indices uniformly at random with replacement from distribution
% whose probability is proportional to the weight vector.
weights = weights/sum(weights);

idxs = zeros(N,1);
cum_weights = cumsum(weights);
for i=1:N
    sampled_ind = find(cum_weights>rand);
    idxs(i) = sampled_ind(1);
end