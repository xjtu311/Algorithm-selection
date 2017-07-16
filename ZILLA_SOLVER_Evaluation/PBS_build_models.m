%% compute score first
yDataSCORE=[];
yDataMETRIC=[];
validationTime=[];

for kk = 1:result.NumSolver
    scorefoo = time2score(yDataOrig(:,kk), result.Cutoff, result.ScoreType);
    yDataSCORE=[yDataSCORE,scorefoo];
    validationTime=[validationTime, yDataOrig(:,kk)+xTimeOrig(:,end)+sum(result.PresolverTime(1:2))];
end
yDataSCORE= yDataSCORE + (rand(size(yDataSCORE))-0.5)*1.0e-05;

for kk = 1:result.NumSolver
    metricfoo = time2score(yDataOrig(:,kk), result.Cutoff, result.MetricType);
    yDataMETRIC=[yDataMETRIC,metricfoo];
    
end

result.Transformation =normalData(xDataOrig);
for jj = 1:result.NumSolver
    xDataFinal=formatnormalData(xDataOrig, result.Transformation);
end

%% all predictins, we will compute RMSE later
tmpresult.TestPred=[];
tmpresult.TrainPred=[];
modelIDSAll=[];

%% have cross validation
startIdx = 1;
for l=1:result.NumCrossValidation
    fprintf(strcat(['Cross-validation ', num2str(l), '/', num2str(result.NumCrossValidation), '**************************************\n']));
    endIdx = ceil(l*result.NumInstance/result.NumCrossValidation);
    testIdx = startIdx:endIdx;
    trainIdx = setdiff(1:result.NumInstance, startIdx:endIdx);
    startIdx = endIdx+1;
    
    %% backup solver
    if result.FeatureCutoff >0
        numInst = length(trainIdx);
        for num_tree=1:99
            X = SF(trainIdx,:);
            y = featlabel(trainIdx);
            feattrees{num_tree} = classregtree(X ,y, 'method', 'classification', 'nvartosample', result.FeatNvartosample);
        end
        predTest=[];
        predTrain=[];
        for num_tree=1:99
            predTest(:, num_tree) = str2num(char(eval(feattrees{num_tree}, SF(testIdx,:))));
            predTrain(:, num_tree) = str2num(char(eval(feattrees{num_tree}, SF(trainIdx,:))));
        end
        featureTimeoutTest = mean(predTest,2);
        featureTimeoutTest(featureTimeoutTest>0)=1;% 1 means feature timeout
        featureTimeoutTest(featureTimeoutTest<0)=-1;
        featureTimeoutTrain = mean(predTrain,2);
        featureTimeoutTrain(featureTimeoutTrain>0)=1;% 1 means feature timeout
        featureTimeoutTrain(featureTimeoutTrain<0)=-1;
    else
        featureTimeoutTest = ones(length(testIdx),1)*-1;
        featureTimeoutTrain = ones(length(trainIdx),1)*-1;
    end
    runtimeTable(testIdx(featureTimeoutTest==1),1)=yDataOrig(testIdx(featureTimeoutTest==1), result.Backupsolver);
    solvedTable(testIdx(featureTimeoutTest==1),[2:end])=0;
    
    %% output classification performance
    %     trainClassification=length(find((featureTimeoutTrain-featlabel(trainIdx))==-2))/length(find(featlabel(trainIdx)==1))
    %     testClassification=length(find((featureTimeoutTest-featlabel(testIdx))==-2))/length(find(featlabel(testIdx)==1))
    %         trainClassification=length(find((featureTimeoutTrain-featlabel(trainIdx))==2))/length(find(featlabel(trainIdx)==-1))
    %     testClassification=length(find((featureTimeoutTest-featlabel(testIdx))==2))/length(find(featlabel(testIdx)==-1))
    
    %% first do presolving. Now we need to have cross validation
    modelIDS=trainIdx(yDataOrig(trainIdx,  result.Presolver(1))> result.PresolverTime(1) & yDataOrig(trainIdx, result.Presolver(2))>result.PresolverTime(2) & featureTimeoutTrain ==-1);
    modelIDSAll=[modelIDSAll, testIdx(yDataOrig(testIdx, result.Presolver(1))>result.PresolverTime(1) & yDataOrig(testIdx, result.Presolver(2))>result.PresolverTime(2) & featureTimeoutTest==-1)];
    %pre1
    solvedbypre1=testIdx(find(yDataOrig(testIdx, result.Presolver(1))<=result.PresolverTime(1)));
    runtimeTable(testIdx,2)=min(yDataOrig(testIdx, result.Presolver(1)), result.PresolverTime(1));
    solvedTable(solvedbypre1,[3:end])=0;
    
    %pre2
    solvedbypre2=testIdx(find(yDataOrig(testIdx, result.Presolver(2))<=result.PresolverTime(2)));
    runtimeTable(testIdx,3)=min(yDataOrig(testIdx, result.Presolver(2)), result.PresolverTime(2));
    solvedTable(solvedbypre2,[4:end])=0;
    
    % feature
    runtimeTable(testIdx,4)=xTime(testIdx,end);
    tmpTrainPred=[];
    tmpTestPred=[];
    %totally. 2) treat them as 0 in the normalized data.
    
    
    for jj = 1:result.NumSolver
        for kk = (jj+1):result.NumSolver
            trainfile=sprintf('%strain%f', svmpath, rand());
            testfile=sprintf('%stest%f', svmpath, rand());
            modelfile=sprintf('%smodel%f', svmpath, rand());
            resultfile=sprintf('%sresult%f', svmpath, rand());
            % first we need cost data and may be want somehow normalize the
            % cost data and we need label the data as -1 and 1
            % what is the cost
            testcostfoo=yDataSCORE(testIdx,jj)-yDataSCORE(testIdx,kk);
            testlabelfoo=ones(length(testcostfoo), 1);
            testlabelfoo(find(testcostfoo<=0))=-1;
            traincostfoo=yDataSCORE(modelIDS,jj)-yDataSCORE(modelIDS,kk);
            trainlabelfoo=ones(length(traincostfoo), 1);
            trainlabelfoo(find(traincostfoo<=0))=-1;
            
            traincostfoo=CostTransformation(traincostfoo, result.Cost);
            testcostfoo=CostTransformation(testcostfoo, result.Cost);
            
            fprintf(strcat(['Cross-validation ', num2str(l), '/', num2str(result.NumCrossValidation), '**************************************\n']));
            
            if strcmp(result.Predictor, 'SVM')
                % we also need output all the features into a file (with -512 or with -512 as 0)
                fout = fopen(trainfile, 'w');
                for uu=1:length(modelIDS)
                    fprintf(fout, '%d cost:%f ', trainlabelfoo(uu), abs(traincostfoo(uu)));
                    for ww=1:size(xDataFinal,2)
                        fprintf(fout, '%d:%f ', ww, xDataFinal(modelIDS(uu), ww));
                    end
                    fprintf(fout, '\n');
                end
                fclose(fout);
                
                fout = fopen(testfile, 'w');
                for uu=1:length(testIdx)
                    fprintf(fout, '%d cost:%f ', testlabelfoo(uu),  abs(testcostfoo(uu)));
                    for ww=1:size(xDataFinal,2)
                        fprintf(fout, '%d:%f ', ww, xDataFinal(testIdx(uu), ww));
                    end
                    fprintf(fout, '\n');
                    
                end
                fclose(fout);
                
                % we need run svm using input file, model file and predition file
                
                mycmd=sprintf('!%ssvm_learn %s %s %s', svmpath, trainfile, modelfile);
                eval(mycmd);
                
                mycmd=sprintf('!%ssvm_classify %s %s %s', svmpath, trainfile, modelfile, resultfile);
                eval(mycmd);
                predfoo1=csvread(resultfile);
                
                mycmd=sprintf('!%ssvm_classify %s %s %s', svmpath, testfile, modelfile, resultfile);
                eval(mycmd);
                predfoo2=csvread(resultfile);
                
                %% remove tmp file as we do that
                delete(trainfile);
                delete(testfile);
                delete(modelfile);
                delete(resultfile);
            end
            
            if strcmp(result.Predictor, 'DF')
                numInst = length(modelIDS);
                if ~isempty(PBS_option.Nvartosample) && PBS_option.Nvartosample>0
                    result.Nvartosample = PBS_option.Nvartosample;
                else
                    result.Nvartosample = ceil(log2(size(xDataFinal,2))+1); % default for classification forest
                end
                for num_tree=1:result.NumTree
                    r = ceil(numInst.*rand(numInst,1));
                    X = xDataFinal(modelIDS,:);
                    y = trainlabelfoo;
                    cost = abs(traincostfoo);
                    if sum(cost) == 0
                        cost = rand(length(cost),1);
                        y = rand(length(cost),1);
                        y(y>0.5)=1;
                        y(y<=0.5)=-1;
                    end
                    trees{num_tree} = classregtree(X ,y, 'method', 'classification', 'weights', cost, 'nvartosample', result.Nvartosample);
                end
                pred=[];
                pred_test=[];
                for num_tree=1:result.NumTree
                    pred(:, num_tree) = str2num(char(eval(trees{num_tree}, xDataFinal(modelIDS,:))));
                    pred_test(:, num_tree) = str2num(char(eval(trees{num_tree}, xDataFinal(testIdx,:))));
                end
                
                predfoo1 = mean(pred,2);
                predfoo2 = mean(pred_test,2);
                predfoo1(find(predfoo1>0))=1;
                predfoo1(find(predfoo1<0))=-1;
                
                
                predfoo2(find(predfoo2>0))=1;
                predfoo2(find(predfoo2<0))=-1;
            end
            
            if strcmp(result.Predictor, 'DT')
                t = classregtree(xDataFinal(modelIDS,:) ,trainlabelfoo, 'method', 'classification', 'weights', abs(traincostfoo));
                predclass1=eval(t, xDataFinal(modelIDS,:));
                predclass2=eval(t, xDataFinal(testIdx,:));
                predfoo1=str2num(char(predclass1));
                predfoo2=str2num(char(predclass2));
            end
            
            tmpTrainPred=[tmpTrainPred, predfoo1];
            tmpTestPred=[tmpTestPred, predfoo2];
        end
    end
    tmpresult.TrainPred=[tmpresult.TrainPred; tmpTrainPred];
    tmpresult.TestPred=[tmpresult.TestPred; tmpTestPred];
    
end
%% let's see the crossvalidation result
indfoo=1;
for ii = 1: result.NumSolver
    for jj = (ii+1): result.NumSolver
        loss=compLossSVM(yDataMETRIC(:,ii), yDataMETRIC(:,jj), tmpresult.TestPred(:,indfoo));
        fprintf('Relative loss for classification between %d and %d is: %f\n', ii, jj, loss);
        result.PredictItem1(indfoo)= ii;
        result.PredictItem2(indfoo) = jj;
        result.RelativeLoss(indfoo) = loss;
        indfoo=indfoo+1;
    end
end

%% With all the predictions, we want know which solver we need to pick
%% based on validation performance
% the input will be predictions of which one is better n*(n-1)/2 of them
% output the solver should be use
if result.NumSolver>1
    tmpresult.ordered = solverSubsetSelectionSVM(tmpresult.TestPred(modelIDSAll,:), validationTime(modelIDSAll,:), result.Cutoff, result.MetricType);
    PBStake=find(tmpresult.ordered.configure(1,:)==1);
else
    PBStake = [1];
end


tmpresult.TestModel=compRUNTIMESVM(PBStake, tmpresult.TestPred,  yDataOrig);

runtimeTable(:,4)=tmpresult.TestModel;
tmpresult.TEST=sum(runtimeTable .* solvedTable, 2);

[PBSscore, PBStime,  PBSsolved] =time2scoretest(tmpresult.TEST,  yDataOrig, result.Cutoff, result.MetricType);

%%=========================================================================
% now, with the selected sub set of solvers, we can build a real model and
% also make the prediction using all the training data

if result.FeatureCutoff >0
    for num_tree=1:99
        X = SF;
        y = featlabel;
        feattrees{num_tree} = classregtree(X ,y, 'method', 'classification', 'nvartosample', result.FeatNvartosample);
    end
    pred=[];
    for num_tree=1:99
        pred(:, num_tree) = str2num(char(eval(feattrees{num_tree}, SF)));
    end
    predfoo1 = mean(pred,2);
    predfoo1(predfoo1>0)=1;
    predfoo1(predfoo1<0)=-1;
    onemodel.trees=feattrees;
    onemodel.type='DF';
    result.featureTimeout=y;
    result.predictfeatureTimeout=predfoo1;
    PBS.FeatModel=onemodel;
else
    PBS.FeatModel=[];
    result.featureTimeout=realFeatTime*0-1;
    result.predictfeatureTimeout=realFeatTime*0-1;
end

modelIDS=find(yDataOrig(:,result.Presolver(1))> result.PresolverTime(1) & yDataOrig(:,result.Presolver(2))>result.PresolverTime(2) & result.predictfeatureTimeout==-1);

modelIndex=1;
if result.NumSolver>1
    PBS.PredictItem1=result.PredictItem1;
    PBS.PredictItem2=result.PredictItem2;
else
    PBS.PredictItem1=[];
    PBS.PredictItem2=[];
end
for jj = 1:result.NumSolver
    for kk = (jj+1):result.NumSolver
        trainfile=sprintf('%strain_%d_d_%f', svmpath, ii, jj, rand());
        modelfile=sprintf('%smodel_%d_%d_%f', svmpath, ii, jj, rand());
        % first we need cost data and may be want somehow normalize the
        % cost data and we need label the data as -1 and 1
        % what is the cost
        traincostfoo=yDataSCORE(modelIDS,jj)-yDataSCORE(modelIDS,kk);
        trainlabelfoo=ones(length(traincostfoo), 1);
        trainlabelfoo(find(traincostfoo<=0))=-1;
        traincostfoo=CostTransformation(traincostfoo, result.Cost);
        
        if strcmp(result.Predictor, 'SVM')
            % we also need output all the features into a file (with -512 or with -512 as 0)
            fout = fopen(trainfile, 'w');
            for uu=1:length(modelIDS)
                fprintf(fout, '%d cost:%f ', trainlabelfoo(uu), abs(traincostfoo(uu)));
                for ww=1:size(xDataFinal,2)
                    fprintf(fout, '%d:%f ', ww, xDataFinal(modelIDS(uu), ww));
                end
                fprintf(fout, '\n');
            end
            fclose(fout);
            
            % we need run svm using input file, model file and predition file
            
            mycmd=sprintf('!%ssvm_learn %s %s %s', svmpath, trainfile, modelfile);
            eval(mycmd);
            
            %% remove tmp file as we do that
            delete(trainfile);
            onemodel.modelfile=modelfile;
            onemodel.type='SVM';
            onemodel.transformation=result.Transformation;
            PBS.Model{modelIndex}=onemodel;
            
        end
        
        if strcmp(result.Predictor, 'DF')
            numInst = length(modelIDS);
            for num_tree=1:result.NumTree
                r = ceil(numInst.*rand(numInst,1));
                X = xDataFinal(modelIDS(r),:);
                y = trainlabelfoo(r);
                cost = abs(traincostfoo(r));
                trees{num_tree} = classregtree(X ,y, 'method', 'classification', 'weights', cost, 'nvartosample', result.Nvartosample);
            end
            onemodel.trees=trees;
            onemodel.type='DF';
            onemodel.transformation=result.Transformation;
            PBS.Model{modelIndex}=onemodel;
        end
        if strcmp(result.Predictor, 'DT')
            t = classregtree(xDataFinal(modelIDS,:) ,trainlabelfoo, 'method', 'classification', 'weights', abs(traincostfoo));
            onemodel.tree=t;
            onemodel.type='DT';
            onemodel.transformation=result.Transformation;
            PBS.Model{modelIndex}=onemodel;
        end
        modelIndex=modelIndex+1;
    end
end
