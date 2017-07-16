function [transformation, xOutput, yOutput, catOutput, catDomainsOutput] = buildCleanData(model_options, namesX, xTrain, yTrain, censTrain, xValid, yValid, censValid, cat, catDomains, numPcaComponents, linearSize, doQuadratic, maxModelSize, responseTransformation, useValid)
% This takes as inputs the raw data xTrain, xValid, yTrain, yValid
% available at training time as well as the names of X and the number
% of parameters. Also some parameters for feature selection.
% It then does the following:
% 0. Move all categorical columns of xTrain and xValid to come first.
% 1. Determine constant columns in xTrain,
%    remove them from xTrain, xValid, and namesX.
%    Decrease numparams if one or more of them are constants.
% 2. Split Xs and names into parameters and features
% 3. Call Lin's feature selection code to select the most informative
%    features (no params involved yet)
% 4. Combine these important features with the parameters
% 5. Normalize combinations w.r.t. the training data
% 6. Again, call feature selection software by KLB to select the
%    most informative combined features
%    (this 2-way process is necessary since memory otherwise explodes)
% 7. Transformation of the response variable Y, e.g. log10
% 8. Return now hopefully clean features and response.
warning off;
if isempty(xValid)
    xValid = xTrain;
    yValid = yTrain;
    censValid = censTrain;
end

%% === Move the categorical columns first in xTrain and xValid.
transformation.cat = cat;
cont = setdiff(1:size(xTrain,2), cat);

tmpCat = xTrain(:,cat);
tmpCont = xTrain(:,cont);
xTrain = [tmpCat, tmpCont];
transformation.cat_idx_mapping = 1:length(cat); % semantic: cat(cat_idx_mapping) = finalCat
tmpNamesCat = namesX(cat);
tmpNamesCont = namesX(cont);
namesX = [tmpNamesCat; tmpNamesCont];

tmpCat = xValid(:,cat);
tmpCont = xValid(:,cont);
xValid = [tmpCat, tmpCont];

numCatParams = length(cat);

if useValid
    xTrainAll=[xTrain;xValid];
else
    xTrainAll=xTrain;
end

%% === 1. Transformation of the response variable Y, e.g. log10
transformation.doQuadratic = doQuadratic;

yTrain = transformResponse(responseTransformation, yTrain);
yValid = transformResponse(responseTransformation, yValid);
transformation.responseTransformation = responseTransformation;

%% 2. Determine constant columns in xTrain (also including xValid),
%    remove them from xTrain, xValid, xTest, and namesX.
%    Decrease numparams if one or more of them are constants.
constants = determine_constants(xTrainAll, 1);
kept_columns = setdiff(1:size(xTrain,2), constants);
xTrain = xTrain(:,kept_columns);
xValid = xValid(:,kept_columns);
xTrainAll = xTrainAll(:,kept_columns);
namesX = namesX(kept_columns);
transformation.cat_idx_mapping = transformation.cat_idx_mapping(kept_columns(find(kept_columns <= numCatParams))); %#ok<FNDSB>
numCatParams= length(transformation.cat_idx_mapping);
transformation.kept_columns = kept_columns;
transformation.numCatParamsAfterRemovingConstants = numCatParams;

%% 3. Split Xs and names into categorical and numerical
% we treat numerical parameters the same as other features
paramsCatTrain = xTrain(:,1:numCatParams);
xTrain = xTrain(:,numCatParams+1:end);
xTrainAll = xTrainAll(:,numCatParams+1:end);
paramsCatValid = xValid(:,1:numCatParams);
xValid = xValid(:,numCatParams+1:end);
transformation.namesX=namesX;
% namesCatParams = namesX(1:numCatParams);
% namesX = namesX(numCatParams+1:end);


%% 4. PCA (if we do PAC, we can not deal with default values)
if numPcaComponents>0
    fprintf('We do not support that! It is in todo list\n');
    fprintf('Becasue we want introduce product of parameters and features \n');
    %       numActualPCAComponents = min([numPcaComponents, size(xTrainAll, 2)]);
    %       [pccoeff, pcVec] = pca(xTrainAll,numActualPCAComponents);
    %       xTrain = xTrain*pcVec;
    %     xValid = xValid*pcVec;
    %     xTrainAll = xTrainAll*pcVec;
    %     transformation.pcVec = pcVec;
else
    transformation.pcVec = eye(size(xTrain,2));
end

%% 5. Linear forward selection to select the most informative features, using normalized data.

% Normalize the data, handling broken features.
% If an entry is one of the predefined problematic values (-512 or -1024),
% then normalization will not count those entries, and the entries will
% be set to zero after normalization.
xTrainDefault=xTrain*0+1;
xTrainDefault(find(xTrain==-512 | xTrain==-1024))=2;
xValidDefault=xValid*0+1;
xValidDefault(find(xValid==-512 | xValid==-1024))=2;
xTrainAllDefault=xTrainAll*0+1;
xTrainAllDefault(find(xTrainAll==-512 | xTrainAll==-1024))=2;

if size(xTrain,2)>linearSize % otherwise, just use all columns
    [nonconstant, scale, bias] = determine_transformation(xTrainAll, xTrainAllDefault, 1);
    if length(nonconstant) < size(xTrain,2)
        error 'At this point, there cannot be any more constant columns; we just removed them above.'
    end
    
    % Make a temporary, normalized, copy of features for linear forward selection
    tmpXTrain = (xTrain - repmat(bias, [size(xTrain,1),1])) ./ repmat(scale, [size(xTrain,1),1]);
    tmpXValid = (xValid - repmat(bias, [size(xValid,1),1])) ./ repmat(scale, [size(xValid,1),1]);
    
    % Set broken features to zero.
    tmpXTrain(find(xTrainDefault> 1))=0;
    tmpXValid(find(xValidDefault>1))=0;
    
    % Add categorical parameters.
    tmpXTrain=[paramsCatTrain, tmpXTrain];
    tmpXValid=[paramsCatValid, tmpXValid];
    
    % Change cat and catDomains to reflect omission of constants and fwd selection.
    updatedCat = 1:numCatParams;
    updatedCatDomains = catDomains(transformation.cat_idx_mapping);
    %=== Assertion to ensure domains are correct.
    for i=1:length(updatedCat)
        for j=1:size(tmpXTrain,1)
            assert(check_member(tmpXTrain(j,i), updatedCatDomains{i}));
        end
    end
    
    % Now, do forward selection!
    [setOfFeatureIndices, bestIndex] = getImportantFeatures(model_options, tmpXTrain, yTrain, censTrain, tmpXValid, yValid, censValid, updatedCat, updatedCatDomains, namesX, linearSize);
    feature_indices = setOfFeatureIndices{bestIndex};
    
    xTrain = xTrain(:,feature_indices(find(feature_indices>numCatParams))-numCatParams);
    xTrainAll = xTrainAll(:,feature_indices(find(feature_indices>numCatParams))-numCatParams);
    xValid = xValid(:,feature_indices(find(feature_indices>numCatParams))-numCatParams);
    
    xTrainDefault=xTrainDefault(:, feature_indices(find(feature_indices>numCatParams))-numCatParams);
    xTrainAllDefault=xTrainAllDefault(:, feature_indices(find(feature_indices>numCatParams))-numCatParams);
    xValidDefault=xValidDefault(:, feature_indices(find(feature_indices>numCatParams))-numCatParams);
    
    transformation.linearFeatureIndices = feature_indices(find(feature_indices>numCatParams));
    transformation.linearFeatureNames = namesX(feature_indices(find(feature_indices>numCatParams)));
    transformation.linearCatParamsNames = namesX(feature_indices(find(feature_indices<=numCatParams)));
    transformation.linearCatParamsIndices=feature_indices(find(feature_indices<=numCatParams));
    transformation.cat_idx_mapping = transformation.cat_idx_mapping(feature_indices(find(feature_indices<=numCatParams)));
    transformation.linearNamesX=namesX(feature_indices);
    transformation.linearNumCatParams=length(find(feature_indices <= numCatParams));
else
    transformation.linearFeatureIndices =  numCatParams+1:length(namesX);
    transformation.linearCatParamsIndices= 1:numCatParams;
    transformation.linearNamesX=namesX;
    transformation.linearNumCatParams=numCatParams;
    transformation.linearFeatureNames = namesX(transformation.linearFeatureIndices);
    transformation.linearCatParamsNames = namesX(transformation.linearCatParamsIndices);
end
numCatParams=length(transformation.cat_idx_mapping);

%% 6. Combine these important features with the parameters
%  TODO: Build clean interface allowing to specify the desired order of the
%  parameters and the instance features. Currently, this is a HACK.
xTrain = buildBasisFunctions(xTrain, transformation.linearFeatureNames, doQuadratic);
xTrainDefault = buildBasisFunctions(xTrainDefault, transformation.linearFeatureNames, doQuadratic);
xValid = buildBasisFunctions(xValid, transformation.linearFeatureNames, doQuadratic);
xValidDefault = buildBasisFunctions(xValidDefault, transformation.linearFeatureNames, doQuadratic);
xTrainAll = buildBasisFunctions(xTrainAll, transformation.linearFeatureNames, doQuadratic);
[xTrainAllDefault, expandedNamesX] = buildBasisFunctions(xTrainAllDefault, transformation.linearFeatureNames, doQuadratic);

%% 7. Normalize combinations w.r.t. the training data

[nonconstant, scale, bias] = determine_transformation(xTrainAll, xTrainAllDefault, 1);
xTrain = xTrain(:,nonconstant);
xValid = xValid(:,nonconstant);
%xTrainAll = xTrainAll(:,nonconstant); % xTrainAll not used afterwards
xTrainDefault = xTrainDefault(:,nonconstant);
xValidDefault = xValidDefault(:,nonconstant);
%xTrainAllDefault = xTrainAllDefault(:,nonconstant); % xTrainAllDefault not used afterwards
xTrain = (xTrain - repmat(bias, [size(xTrain,1),1])) ./ repmat(scale, [size(xTrain,1),1]);
xValid = (xValid - repmat(bias, [size(xValid,1),1])) ./ repmat(scale, [size(xValid,1),1]);
transformation.nonconstant = nonconstant;
transformation.scale=scale;
transformation.bias=bias;
expandedNamesX = expandedNamesX(nonconstant);

%% 8. Again, call forward feature selection (on the expanded features) to select the most informative combined features
%    (this 2-way process is necessary since memory otherwise explodes)

xTrain(find(xTrainDefault> 1))=0;
xValid(find(xValidDefault>1))=0;
xTrain=[paramsCatTrain(:,transformation.linearCatParamsIndices), xTrain];
xValid=[paramsCatValid(:,transformation.linearCatParamsIndices), xValid];
allNamesX = [transformation.linearCatParamsNames;expandedNamesX];

if (maxModelSize >0 && maxModelSize<size(xTrain,2))
    % Change cat and catDomains to reflect omission of constants and fwd selection.
    updatedCat = 1:numCatParams;
    updatedCatDomains = catDomains(transformation.cat_idx_mapping);
    %=== Assertion to ensure domains are correct.
    for i=1:length(updatedCat)
        for j=1:size(xTrain,1)
            assert(check_member(xTrain(j,i), updatedCatDomains{i}));
        end
    end
    
    [setOfFeatureIndices, bestIndex] = getImportantFeatures(model_options, xTrain, yTrain, censTrain, xValid, yValid, censValid, updatedCat, updatedCatDomains, allNamesX, maxModelSize);
else
    setOfFeatureIndices={[1:size(xTrain,2)]};
    bestIndex=1;
end

transformation.setOfFinalFeatureIndices = setOfFeatureIndices;
transformation.bestIndexForFinalFeatures = bestIndex;
feature_indices = setOfFeatureIndices{bestIndex};
transformation.finalFeatureIndices=feature_indices(find(feature_indices>numCatParams));
transformation.finalCatParamsIndices=feature_indices(find(feature_indices<=numCatParams));
transformation.finalFeatureNames=allNamesX(feature_indices(find(feature_indices>numCatParams)));
transformation.finalCatParamsNames=allNamesX(feature_indices(find(feature_indices<=numCatParams)));
transformation.cat_idx_mapping = transformation.cat_idx_mapping(feature_indices(find(feature_indices<=numCatParams)));
transformation.finalNamesX=[transformation.finalCatParamsNames; transformation.finalFeatureNames];

tmpCatxTrain = xTrain(:,transformation.finalCatParamsIndices);
tmpxTrain = xTrain(:,transformation.finalFeatureIndices);
xOutput=[tmpCatxTrain, tmpxTrain];

yOutput = yTrain;
catOutput = 1:length(find(feature_indices<=numCatParams));
catDomainsOutput = catDomains(transformation.cat_idx_mapping);

%=== Assertion to ensure domains are correct.
for i=1:length(catDomainsOutput)
    for j=1:size(xOutput,1)
        assert(check_member(xOutput(j,i), catDomainsOutput{i}));
    end
end


function result = check_member(value, cell_array)
result = (value <= length(cell_array));
% for i=1:length(cell_array)
%     if strcmp(cell_array{i}, num2str(value))
%         result = true;
%         return;
%     end
% end
% result = false;