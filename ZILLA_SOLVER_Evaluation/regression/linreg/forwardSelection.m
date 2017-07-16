function subs = forwardSelection(xTrain, yTrain, xValid, yValid, allnames, maxSize)
%function subs = subsetSelection(df, maxSize, type)
%
%Do subset selection on df, using the forward selection heuristic, looking
%at subsets of increasing size.
%Returns a struct array of size 1..maxSize, containing selected features
%vector (colidx), rmse and mae of the model (on training data)


%display (sprintf('Forward selection started: %d features total',size(df.traindata.features,2)));
X = ones(size(xTrain,1),1);
delta=1e-2; % = 1e-16;
A = X'*X+delta;
Ainv = inv(A);
XtY = X' * yTrain;
compute_training_RMSE = false;

%test = feval(df.backtransform, df.validdata.response, df);
test = yValid;

features = [];
num_train_instances = size(xTrain,1);
num_test_instances = size(xValid,1);
num_features = size(xTrain,2);
maxIteration = min(maxSize,num_features);
Xtest = [ones(num_test_instances,1) zeros(num_test_instances,maxIteration+1)];

% for every set size
for t=1:maxIteration
    
    minRMSE = inf;
    best_subset = 0;
    
    if (mod(t,5)==0)
        % display (sprintf('  forward selection: subset size %d of %d', [t maxIteration]));
    end
    
    % for each feature not in the set
    for tt=1:size(xTrain,2) %- size(features,2)
        if (~isempty(find(features==tt)))
            continue;
        end
        
        features(t)=tt;
        
        % construct the model with this feature added
        model = modelAddFeature(X, xTrain(:,tt), yTrain, A, Ainv, delta, XtY);
        
        % record the RMSE and MAE on validation data for this feature
        Xtest(:,t+1) = xValid(:,tt);
        preds = Xtest*[model; zeros(maxIteration - t+1,1)];
        
        % are we sure we want to do this?  This means we're doing model
        % selection based on a different error measure from the one we
        % optimize in the regression...
        %		preds = feval(df.backtransform, preds, df);
        
        result(tt).MAE = mean(abs(preds-test));
        result(tt).RMSE = sqrt(mean( (preds-test).^2 ));
        
        if (compute_training_RMSE)
            preds = [X xTrain(:,tt)]*model;
            %preds = feval(df.backtransform, preds, df);
            %result(tt).trainingRMSE = sqrt(mean((preds-feval(df.backtransform, df.traindata.response, df)).^2));
            result(tt).trainingRMSE = sqrt(mean((preds-yTrain).^2));
        end
        
        result(tt).features = features;
        result(tt).method = 'forwardSelection';
        minRMSE = min(minRMSE, result(tt).RMSE);
        if (minRMSE == result(tt).RMSE)
            best_subset = tt;
        end
    end
    
    % update the inverse; global variables
    %% this line is modified by lin at 2/4/2006 and put some lines below
%     [model A Ainv XtY] = modelAddFeature(X, xTrain(:,best_subset), yTrain, A, Ainv, delta, XtY);
%     
    %% keep this two line always
    X = [X xTrain(:,best_subset)];
    Xtest(:,t+1) = xValid(:,best_subset);
    
    %%% replace those line with 3 lines up
        A = X' * X + delta * eye(size(X,2));
        Ainv=pinv(A);
        XtY = X' * yTrain;
        model=Ainv * XtY;
        preds = Xtest*[model; zeros(maxIteration - t+1,1)];
        result(best_subset).MAE = mean(abs(preds-test));
        result(best_subset).RMSE = sqrt(mean( (preds-test).^2 ));
    
    % pick the feature with the best RMSE
    subs(t) = result(best_subset);
    clear result;
    features = subs(t).features;
    
end

% display 'Forward selection finished';