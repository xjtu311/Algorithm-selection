function seeds = getTrainSeedsForInstanceNumbers(func, Theta_idx, instanceNumbers, increaseCount)
% Pick every seed uniformly at random from the user-provided (or generated)
% seeds for the instance. Blocking introduces bias in the following sense:
% when all points in the initial LHD use the same seed, we learn really
% well about that seed, but not about the behaviour with all seeds.
% Think about blocking on instances: we wouldn't want to just use 1
% instance for the whole initial LHD, but instead learn about performance
% with all instances. So not even blocking on instances seems like a good
% idea here...
% This is similar in intuition to the result that
% picking "N instances, 1 run each" leads to a minimal expected
% variance estimator for the mean.
for i=1:length(instanceNumbers)
    seeds(i) = func.seeds(instanceNumbers(i), ceil(size(func.seeds,2)*rand));
end
