function [setOfFeatureIndices, bestIndex] = feature_select(xTrain, yTrain, xValid, yValid, allnames, maxSize)
%=== Currently only calls forward selection and returns all its subsets,
%=== plus the index for which where validation set RMSE was lowest.

subs = forwardSelection(xTrain, yTrain, xValid, yValid, allnames, maxSize);
rmses = [];
for i=1:length(subs)
    setOfFeatureIndices{i} = subs(i).features;
    rmses(i) = subs(i).RMSE;
end

[tmp, bestIndex] = min(rmses);
% figure;
% plot([1:length(subs)], rmses, 'r.');
% fprintf('The minimal rmse with %d features is %f \n', bestIndex, tmp);