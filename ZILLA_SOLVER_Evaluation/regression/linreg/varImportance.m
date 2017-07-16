function vimp1=varImportance(xTrain, yTrain, xValid, yValid, tansformation, modelType)
%function vimp=varImportance(df)
%
%Do variable importance, by collecting stats for each leave-one variable
%model.
n=size(xTrain,2);

fullmodel = learnModel(xTrain, yTrain,  [], [],  modelType, []);
[yPredValid, yPredVarValid] = applyModel(fullmodel, xValid);
fullRMSE = sqrt(mean((yPredValid-yValid).^2 ));

for i=1:n
    colidx=((1:n)~=i);
    tmpxTrain = xTrain(:,colidx);
    tmpxValid = xValid(:,colidx);
    tmpmodel = learnModel(tmpxTrain, yTrain,  [], [],  modelType, []);
    [tmpyPredValid, tmpyPredVarValid] = applyModel(tmpmodel, tmpxValid);
    tmpRMSE = sqrt(mean((tmpyPredValid-yValid).^2 ));
    imp(i) = tmpRMSE - fullRMSE;
end

%normalize
imp = 100 * (imp ./ max(imp));
[tmp idx] = sort(imp);
% build the variable importance structure
vimp = struct('feature', tansformation.finalNamesX(idx)','importance', num2cell(imp(idx)));
% vimp1.Aid=tansformation.Aid(idx);
% vimp1.Bid=tansformation.Bid(idx);
vimp1.idx=idx;
vimp1.impotance=imp(idx);

figure;
barh([vimp.importance]);
set(gca,'YTickLabel',char({vimp.feature}'), 'FontSize', 14);
title (sprintf('Relative Importance of Variables in Subset of Size %d (validation data)',size([vimp.importance],2)), 'FontSize', 14);
xlabel ('Relative Importance', 'FontSize', 14);

