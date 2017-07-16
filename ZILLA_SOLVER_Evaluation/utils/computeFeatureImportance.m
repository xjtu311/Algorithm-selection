function computeFeatureImportance(dataset, it, outputfile)
mypath='./data/';
if isdeployed
    it = str2double(it);
    mypath='/ubc/cs/project/arrow/projects/FORLIN/MipHydra/HYDRA/SATZILLA/data/';
end
% some parameters you may want to chance

nLSSOLVER = 0; %% indu category
pretimeLimit=0;
nSOLVER = it;
cutoff =3600;
responseTransformation = 'identity';
take = [1:nSOLVER];
lstake=[];
% some information about features
reference = take;
goodfeatures=[1:40];
timefeatures=[127:129];
goodfeatures=[1:139];
timefeatures=[140:142];
% goodfeatures=[1:94];
% timefeatures=[140:142];
File.name= sprintf('%s%s-solver.csv', mypath, dataset);
testFile.name= sprintf('%s%s-solver.csv', mypath, dataset);
oFile.name= sprintf('%s%s-default.csv', mypath, dataset);
testoFile.name= sprintf('%s%s-default.csv', mypath, dataset);
featFile.name= sprintf('%s%s-feat.csv', mypath, dataset);
testfeatFile.name= sprintf('%s%s-feat.csv', mypath, dataset);

% ++++++++++++++++++++++++++++++++++++++++
useSCORE = 1;
useSERIAL =0;
% nCLASS = 3; % if we only have one class, then it means only take the first class
serial =[];
numOutputs = 5;
removeCapped = 1;
eps=10^-5;
linearSize = 40;
doQuadratic = 1;
maxModelSize = 20;
maxMESize = 20;
maxLSMESize = 15;
numParams =0;
numPcaComponents = 0; %(set to 0: PCA turned off)
modelType = 'LR';
miniResponse = 0.005;
alltake = [take,lstake];

% step one read data
% make sure we also read the instance id data
[namesY, namesX, yData.SAT, solution, xData.SAT] = readRawDataHydra(File.name, oFile.name, featFile.name, numOutputs, nSOLVER+nLSSOLVER, miniResponse);
[testnamesY, testnamesX, testyData.SAT, testsolution, testxData.SAT] = readRawDataHydra(testFile.name, testoFile.name, testfeatFile.name, numOutputs, nSOLVER+nLSSOLVER, miniResponse);

trainfeattime=max(xData.SAT(:, timefeatures), 0);
testfeattime=max(testxData.SAT(:, timefeatures), 0);
yData.SAT = min(yData.SAT, cutoff); 
xData.SAT=[xData.SAT, sum(trainfeattime,2)];
testyData.SAT = min(testyData.SAT, cutoff); 
testxData.SAT=[testxData.SAT,sum(testfeattime,2)];
% select part of feature set (some of features may be bad constance)
namesX=namesX(goodfeatures);
% always good to have something as back up
xData.ALLORIG=[xData.SAT];
yData.ALLORIG=[yData.SAT];

testxData.ALLORIG=[testxData.SAT];
testyData.ALLORIG=[testyData.SAT];
% we only care about test performance (in training phase, training file is the same as test file)
% my new structure will be like this. For test data, we have two strucure
% for recording runtime results (PAR score will be computed based on run time).
% one is runtime table, the second is solved table
% realruntime for hydra will be the dot product of these two
runtimeTable=zeros(size(testyData.SAT,1), 4);
solvedTable=runtimeTable+1;

preSolver1 = 2;
preTime1 = 5;
preSolver2 = 3;
preTime2 = 5;

hydralinux;
% write out put into outfile

