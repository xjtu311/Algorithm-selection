function [yPred, yPredVar] = applyModel(model, X)
maxPred = 400000;
minPred = -400000;
switch model.type
    case 'BLR'
        % for Bayesian linear regression
        [yPred, yPredVar] = bayes_fwd(model.mu, model.Sigma, model.beta, X);
        yPred = min(max(yPred, minPred), maxPred);
    case 'LR'
        % for linear regression
        [yPred, yPredVar] = LRT(model, X);
        yPred = min(max(yPred, minPred), maxPred);
        
    case 'matlab-regtree'
        yPredMean = treeval(model.T, X);
        yPredVar = zeros(length(yPredMean),1);
        
    case 'fh-regtree'
        yPredMean = fh_simple_one_treeval(model.T, X);
        if model.options.logModel
            yPredMean = log10(yPredMean);
        end
        yPredVar = zeros(length(yPredMean),1);
        
    case 'rf'
        treemeans = zeros(size(X,1), length(model.module));
        for m=1:length(model.module)
            treemeans(:,m) = fh_simple_one_treeval(model.module{m}.T, X);
        end
        if model.options.logModel
            yPred = mean(log10(treemeans),2);
            yPredVar = var(log10(treemeans),0,2);
        else
            yPred = mean(treemeans,2);
            yPredVar = var(treemeans,0,2);
        end
    otherwise
        error ('No such model type defined yet!');
end
