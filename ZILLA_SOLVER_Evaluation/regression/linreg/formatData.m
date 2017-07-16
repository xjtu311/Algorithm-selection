function [xData, yData, cat, catDomains] = formatData(xData, yData, transformation)

%=== Move the categorical columns first.
cont = setdiff(1:size(xData,2), transformation.cat);
tmpCat = xData(:,transformation.cat);
tmpCont = xData(:,cont);
xData = [tmpCat, tmpCont];

%=== Transform response variable.
yData = transformResponse(transformation.responseTransformation, yData);


%=== Remove constant features.
xData = xData(:,transformation.kept_columns);

%=== Cat.parameters 
paramsCatxData=xData(:,[1:transformation.numCatParamsAfterRemovingConstants]);

%=== Split into parameters and instance features.
xData = xData(:,transformation.numCatParamsAfterRemovingConstants+1:end);


%=== PCA. (no PCA at this point)
% xData = xData*transformation.pcVec;



%=== Feature selection for linear instance features.
foo=[paramsCatxData, xData];
xData = foo(:,transformation.linearFeatureIndices);
paramsCatxData = foo(:, transformation.linearCatParamsIndices);
%=== deal with default feature
xDataDefault=xData*0+1;
xDataDefault(find(xData==-512 | xData==-1024))=2;
%=== Build basis functions.
%=== For data without parameters, this may include building quadratic features etc.
%=== For data with parameters, it also builds combinations of parameters and instance features.
[xData, NamesX, Aid, Bid] = buildBasisFunctions(xData, transformation.linearFeatureNames, transformation.doQuadratic);
[xDataDefault, NamesX, Aid, Bid] = buildBasisFunctions(xDataDefault, transformation.linearFeatureNames, transformation.doQuadratic);

%=== Normalize data.
xData = xData(:, transformation.nonconstant);
xData = (xData - repmat(transformation.bias, [size(xData,1),1])) ./ repmat(transformation.scale, [size(xData,1),1]);
xDataDefault = xDataDefault(:, transformation.nonconstant);
xData(find(xDataDefault>1))=0;

finalxData = [paramsCatxData, xData];
%=== Final feature selection to build small model.
% finalFeatureIndices = transformation.setOfFinalFeatureIndices{transformation.bestIndexForFinalFeatures};
% xData = finalxData(:,finalFeatureIndices);

%=== Put categorical parameters first.
tmpCatxTrain = finalxData(:,transformation.finalCatParamsIndices);
tmpxTrain = finalxData(:,transformation.finalFeatureIndices);
xData=[tmpCatxTrain, tmpxTrain];

t = any(any(isnan(xData),2));
t = t | any(any(isinf(xData),2));
t = t | any(isnan(yData));
t = t | any(isinf(yData));
if any(t)
    error 'Empty values and infinity not allowed in either X, y.'
end