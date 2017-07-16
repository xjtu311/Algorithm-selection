% this code is for select right cutoff point for sat and unsat 
% 1. we train three models for sat unsat and unconditional
% 2. we choose the optimal threshold fro using sat model unsatmodel based
% on the validation data interms of minimze the RMSE
% 3. for the new instances, we select the right model based on the
% classification score threeshold we set uo in previous step.

%========================================================================
% 1. Configure parameters.
%========================================================================
% trainFilename = 'C:\Temp\classify\newsat\kcnfsFIXALLDATA1-train.csv';
% validFilename = 'C:\Temp\classify\newsat\kcnfsFIXALLDATA1-valid.csv';
% testFilename = 'C:\Temp\classify\newsat\kcnfsFIXALLDATA1-test.csv';
% trainFilename = 'C:\Temp\classify\newsat\kcnfsVARALLDATA1-train.csv';
% validFilename = 'C:\Temp\classify\newsat\kcnfsVARALLDATA1-valid.csv';
% testFilename = 'C:\Temp\classify\newsat\kcnfsVARALLDATA1-test.csv';
% trainFilename = 'C:\Temp\classify\newsat\satzooSWALLDATA2-train.csv';
% validFilename = 'C:\Temp\classify\newsat\satzooSWALLDATA2-valid.csv';
% testFilename = 'C:\Temp\classify\newsat\satzooSWALLDATA2-test.csv';
% trainFilename = 'C:\Temp\classify\newsat\oksolverFIXALLDATA1-train.csv';
% validFilename = 'C:\Temp\classify\newsat\oksolverFIXALLDATA1-valid.csv';
% testFilename = 'C:\Temp\classify\newsat\oksolverFIXALLDATA1-test.csv';
% trainFilename = 'C:\Temp\classify\newsat\satzFIXALLDATA1-train.csv';
% validFilename = 'C:\Temp\classify\newsat\satzFIXALLDATA1-valid.csv';
% testFilename = 'C:\Temp\classify\newsat\satzFIXALLDATA1-test.csv';
% trainFilename = 'C:\Temp\classify\newsat\satzooSWALLDATA2-train.csv';
% validFilename = 'C:\Temp\classify\newsat\satzooSWALLDATA2-valid.csv';
% testFilename = 'C:\Temp\classify\newsat\satzooSWALLDATA2-test.csv';
% trainFilename = 'C:\Temp\classify\newsat\sateliteQCPALLDATA2-train.csv';
% validFilename = 'C:\Temp\classify\newsat\sateliteQCPALLDATA2-valid.csv';
% testFilename = 'C:\Temp\classify\newsat\sateliteQCPALLDATA2-test.csv';
% trainFilename = 'C:\Temp\classify\newsat\satzooQCPALLDATA2-train.csv';
% validFilename = 'C:\Temp\classify\newsat\satzooQCPALLDATA2-valid.csv';
% testFilename = 'C:\Temp\classify\newsat\satzooQCPALLDATA2-test.csv';

trainFilename = 'C:\Temp\classify\newsat\SW300-C6-10-7__satelite__1__3600.0-train.csv';
validFilename = 'C:\Temp\classify\newsat\SW300-C6-10-7__satelite__1__3600.0-valid.csv';
testFilename = 'C:\Temp\classify\newsat\SW300-C6-10-7__satelite__1__3600.0-test.csv';
%%% QCP remove instance with nvar less than 100 %%%%%%%%%%
% trainFilename = 'C:\Temp\classify\newsat\sateliteQCPALLDATA3-train.csv';
% validFilename = 'C:\Temp\classify\newsat\sateliteQCPALLDATA3-valid.csv';
% testFilename = 'C:\Temp\classify\newsat\sateliteQCPALLDATA3-test.csv';
% testFilename = validFilename;

numOutputs = 1;
removeCapped = 1;
numParams = 0;
responseTransformation = 'log10';
linearSize = 30;
doQuadratic = 1;
maxModelSize = 80;
numPcaComponents = 0; %(set to 0: PCA turned off)
modelAndTransformationFilename = 'satzooQCPModel';
%modelType = 'BLR';
modelType = 'LR';
miniResponse = 0.01;

%========================================================================
% 2. Read raw training and validation data
%========================================================================
[namesY, namesX, yTrain, cappedTrain, solutionTrain, xTrain] = readRawData(trainFilename, numOutputs, removeCapped, miniResponse);
[namesY, namesX, yValid, cappedValid, solutionValid, xValid] = readRawData(validFilename, numOutputs, removeCapped, miniResponse);


%========================================================================
% 2. classifiy to two classes
%========================================================================
% change class label to stander format
solutionTrain(find(solutionTrain==-1))=2;
solutionValid(find(solutionValid==-1))=2;

% begain iteration (original value was 10)

% train a classifier and test on traindata and validation data
%[w, log_posterior, kernel, kernel_param, X_train] = classifyTrain(xTrain(1:20000,:), solutionTrain(1:20000));
[w, log_posterior, kernel, kernel_param, X_train] = classifyTrain(xTrain, solutionTrain);
[pValid, classValid] = classifyTest(xValid, solutionValid, w,log_posterior, kernel, kernel_param, X_train);
[pTrain, classTrain] = classifyTest(xTrain, solutionTrain, w,log_posterior, kernel, kernel_param, X_train);
cerrorTrain=size(find(solutionTrain~=classTrain),1);
cerrorValid=size(find(solutionValid~=classValid),1);

% seprate original data into two class based on class label
xTrainSat=xTrain(find(solutionTrain(:,1) ==1),:);
yTrainSat=yTrain(find(solutionTrain(:,1)==1));
xTrainUnsat=xTrain(find(solutionTrain(:,1) ==2),:);
yTrainUnsat=yTrain(find(solutionTrain(:,1) ==2));

xValidSat=xValid(find(solutionValid(:,1) ==1),:);
yValidSat=yValid(find(solutionValid(:,1)==1));
xValidUnsat=xValid(find(solutionValid(:,1) ==2),:);
yValidUnsat=yValid(find(solutionValid(:,1) ==2));

%========================================================================
% 3. Transform raw data into nice data for ML and remember the
%    transformation.
%========================================================================
% find the best basis function using train data and validation data
transformationSat = buildCleanData(xTrainSat, xValidSat, yTrainSat, yValidSat, namesX, numParams, numPcaComponents, linearSize, doQuadratic, maxModelSize, responseTransformation, xTrain);
transformationSat.removeCapped=removeCapped;
transformationSat.miniResponse=miniResponse;

transformationUnsat = buildCleanData(xTrainUnsat, xValidUnsat, yTrainUnsat, yValidUnsat, namesX, numParams, numPcaComponents, linearSize, doQuadratic, maxModelSize, responseTransformation, xTrain);
transformationUnsat.removeCapped=removeCapped;
transformationUnsat.miniResponse=miniResponse;

transformation = buildCleanData(xTrain, xValid, yTrain, yValid, namesX, numParams, numPcaComponents, linearSize, doQuadratic, maxModelSize, responseTransformation, xTrain);
transformation.removeCapped=removeCapped;
transformation.miniResponse=miniResponse;

%========================================================================
% 4. Apply transformation to the raw data.
%========================================================================
[xFinal, yFinal] = formatData(xTrain, yTrain, transformation);
[xFinalSat, yFinalSat] = formatData(xTrainSat, yTrainSat, transformationSat);
[xFinalUnsat, yFinalUnsat] = formatData(xTrainUnsat, yTrainUnsat, transformationUnsat);
%========================================================================
% 5. Build and save model.
%========================================================================
modelSat = learnModel(xFinalSat, yFinalSat, modelType);
%save (modelAndTransformationFilename, 'model', 'transformation');
modelUnsat = learnModel(xFinalUnsat, yFinalUnsat, modelType);
model = learnModel(xFinal, yFinal, modelType);

[xValidFinalSat, yValidFinalSat] = formatData(xValid, yValid, transformationSat);
[xValidFinalUnsat, yValidFinalUnsat] = formatData(xValid, yValid, transformationUnsat);
[xValidFinal, yValidFinal] = formatData(xValid, yValid, transformation);
% compute reponse for validation data
[yPredSat, yPredVarSat] = applyModel(modelSat, xValidFinalSat);
[yPredUnsat, yPredVarUnsat] = applyModel(modelUnsat, xValidFinalUnsat);
[yPred, yPredVar] = applyModel(model, xValidFinal);
% find the best cutoff point for using sat model or unsat model
tmprmseSat = [];
for i = 0.001:0.0010:1
    tmpvalidpreds = yPred;
    tmpvalidpreds(find(pValid(:,1) > (1-i))) = yPredSat(find(pValid(:,1) > (1-i)));
    tmprmseSat=[tmprmseSat, sqrt(mean((tmpvalidpreds-yValidFinalSat).^2 ))];
end
satCutoff = max(1-0.001*find(tmprmseSat==min(tmprmseSat)));
tmprmseUnsat = [];
for i = 0.0001:0.00010:1
    tmpvalidpreds = yPred;
    tmpvalidpreds(find(pValid(:,1) < i)) = yPredUnsat(find(pValid(:,1) < i));
    tmprmseUnsat=[tmprmseUnsat, sqrt(mean((tmpvalidpreds-yValidFinalSat).^2 ))];
end
unsatCutoff = min(0.0001*find(tmprmseUnsat==min(tmprmseUnsat)));
%relabel the class label for trainnig data and validation data


%=====================================================
%            Test                                    %
%=====================================================

[namesY, namesX, yTest, cappedTest, solutionTest, xTest] = readRawData(testFilename, numOutputs, transformationSat.removeCapped, transformationSat.miniResponse);
solutionTest(find(solutionTest==-1))=2;
[pTest, classTest] = classifyTest(xTest, solutionTest, w,log_posterior, kernel, kernel_param, X_train);
cerrorTest=size(find(solutionTest~=classTest),1);
%========================================================================
% 2. Apply transformation.
%========================================================================
[xTestFinalSat, yTestFinalSat] = formatData(xTest, yTest, transformationSat);
[xTestFinalUnsat, yTestFinalUnsat] = formatData(xTest, yTest, transformationUnsat);
[xTestFinal, yTestFinal] = formatData(xTest, yTest, transformation);

%========================================================================
% 3. Apply model.
%========================================================================
[yPredTestSat, yPredVarTestSat] = applyModel(modelSat, xTestFinalSat);
[yPredTestUnsat, yPredVarTestUnsat] = applyModel(modelUnsat, xTestFinalUnsat);
[yPredTest, yPredVarTest] = applyModel(model, xTestFinal);
testPreds = yPredTest;
testPreds(find(pTest(:,1) > satCutoff)) = yPredTestSat(find(pTest(:,1) > satCutoff));
testPreds(find(pTest(:,1) < unsatCutoff)) = yPredTestUnsat(find(pTest(:,1) < unsatCutoff));
%%% add minimal runtime bound %%%
% testPreds(find(testPreds < log10(miniResponse)))=log10(miniResponse);
testRMSE= sqrt(mean((testPreds-yTestFinalSat).^2 ));

magicpreds = yPredTestUnsat;
magicpreds(find(solutionTest==1)) = yPredTestSat(find(solutionTest==1));

mixPreds = yPredTestSat .* pTest(:,1) + yPredTestUnsat .* pTest(:,2);
naiveRMSE= compRMSE(mixPreds, yTestFinalSat);

%   magicpreds(find(magicpreds < log10(miniResponse)))=log10(miniResponse);
magicRMSE = sqrt(mean((magicpreds-yTestFinalSat).^2 ));
uncondRMSE = sqrt(mean((yPredTest-yTestFinalSat).^2 ));
satallRMSE = sqrt(mean((yPredTestSat-yTestFinalSat).^2 ));
unsatallRMSE = sqrt(mean((yPredTestUnsat-yTestFinalSat).^2 ));


fprintf('===============================================================\n');
fprintf('                  RESULTS for\n')
fprintf('                  %s\n', trainFilename);
fprintf('===============================================================\n');

fprintf('Train data classififcation error is %f \n', cerrorTrain/size(solutionTrain,1));
fprintf('Validation data classififcation error is %f \n', cerrorValid/size(solutionValid,1));
fprintf('Test data classififcation error is %f \n', cerrorTest/size(solutionTest,1));
fprintf('_______________________________________________________________\n');
fprintf('The sat cutoff point is %f \n', satCutoff);
fprintf('The unsat cutoff point is %f \n', unsatCutoff);
fprintf('Prediction RMSE of mixture model with optimal cutoff %f \n', testRMSE);
fprintf('_______________________________________________________________\n');

fprintf('Prediction RMSE of magic model %f \n', magicRMSE);
fprintf('Prediction RMSE of unconditional model %f \n', uncondRMSE);
fprintf('Prediction RMSE of Naive mixture of expert  %f \n', naiveRMSE);
fprintf('Prediction RMSE of SAT model ONLY %f \n', satallRMSE);
fprintf('Prediction RMSE of UNSAT model only %f \n', unsatallRMSE);
fprintf('===============================================================\n');

mtitle = sprintf('Prediction RMSE of magic model %f \n', magicRMSE);
untitle = sprintf('Prediction RMSE of unconditional model %f \n', uncondRMSE);
ntitle = sprintf('Prediction RMSE of Naive mixture of expert  %f \n', naiveRMSE);
stitle = sprintf('Prediction RMSE of SAT model ONLY %f \n', satallRMSE);
utitle = sprintf('Prediction RMSE of UNSAT model only %f \n', unsatallRMSE);
ttitle = sprintf('Prediction RMSE of mixture model with optimal cutoff %f \n', testRMSE);

plotFig(magicpreds, yTestFinalSat, mtitle);
plotFig(yPredTest, yTestFinalSat, untitle);
plotFig(mixPreds, yTestFinalSat, ntitle);
plotFig(yPredTestSat, yTestFinalSat, stitle);
plotFig(yPredTestUnsat, yTestFinalSat, utitle);
plotFig(testPreds, yTestFinalSat, ttitle);





%     testPreds = yPredTestUnsat;
%     testPreds(find(solutionTest == 1)) = yPredTestSat(find(solutionTest == 1));
%     testRmse= sqrt(mean((testPreds-yTestFinalSat).^2 ));
%     fprintf('Test rmse now is %f \n', testRmse);