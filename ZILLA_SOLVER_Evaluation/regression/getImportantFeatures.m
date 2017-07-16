function [setOfFeatureIndices, bestIndex] = getImportantFeatures(model_options, xTrain, yTrain, censTrain, xValid, yValid, censValid, cat, catDomains, allnames, maxSize)

switch model_options.modelType
    case 'LR'
        subs = forwardSelection(xTrain, yTrain, xValid, yValid, allnames, maxSize);
        for i=1:length(subs)
            setOfFeatureIndices{i} = subs(i).features;
            rmses(i) = subs(i).RMSE;
        end
        [tmp, bestIndex] = min(rmses);
        
    otherwise
        xInTrain = zeros(size(xTrain,1),0);
        xInValid = zeros(size(xValid,1),0);
        chosen = [];
        while size(xInTrain,2) < maxSize
            rmse = inf*ones(size(xTrain,2),1);
            for i=1:size(xTrain,2)
                if ~ismember(i,chosen) 
                    model = learnModel([xInTrain, xTrain(:,i)], yTrain, censTrain, cat, catDomains, 1, model_options, allnames);
                    yPred = applyModel(model, [xInValid, xValid(:,i)], 1, 0, 0);
                    rmse(i) = compRMSE(yPred, yValid);
                end
            end
            [tmp, minidx] = min(rmse);
            chosen = [chosen, minidx];
            xInTrain = [xInTrain, xTrain(:,minidx)];
            xInValid = [xInValid, xValid(:,minidx)];
        end
        setOfFeatureIndices = {chosen};
        bestIndex = 1;
end