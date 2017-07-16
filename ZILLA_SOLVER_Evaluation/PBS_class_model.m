%% compute score first
yDataSCORE=[];
yDataMETRIC=[];
ALLyDataSCORE=[];
ALLyDataMETRIC=[];
validationTime=[];

for kk = 1:result.NumSolver
    scorefoo = time2score(yDataOrig(:,kk), result.Cutoff, result.ScoreType);
    yDataSCORE=[yDataSCORE,scorefoo];
    validationTime=[validationTime, yDataOrig(:,kk)+xTimeOrig(:,end)+sum(result.PresolverTime(1:2))];
    scorefoo = time2score(ALLyData(:,kk), result.Cutoff, result.ScoreType);
    ALLyDataSCORE=[ALLyDataSCORE,scorefoo];
end
yDataSCORE= yDataSCORE + (rand(size(yDataSCORE))-0.5)*1.0e-05;

for kk = 1:result.NumSolver
    metricfoo = time2score(yDataOrig(:,kk), result.Cutoff, result.MetricType);
    yDataMETRIC=[yDataMETRIC,metricfoo];
    metricfoo = time2score(ALLyData(:,kk), result.Cutoff, result.MetricType);
    ALLyDataMETRIC=[ALLyDataMETRIC,metricfoo];
    
end

result.Transformation =normalData(xDataOrig);
PBS.Transformation = result.Transformation;
for jj = 1:result.NumSolver
    xDataFinal=formatnormalData(xDataOrig, result.Transformation);
    ALLxDataFinal=formatnormalData(ALLxData, result.Transformation);
end

%% all predictins, we will compute RMSE later
tmpresult.TestPred=[];
tmpresult.TrainPred=[];
modelIDSAll=[];

%%=========================================================================
% now, with the selected sub set of solvers, we can build a real model and
% also make the prediction using all the training data

if result.FeatureCutoff >0
    for num_tree=1:99
        feattrees{num_tree} = classregtree(xsData, featLabel, 'method', 'classification', 'nvartosample', result.FeatNvartosample);
    end
    pred=[];
    for num_tree=1:99
        pred(:, num_tree) = str2num(char(eval(feattrees{num_tree}, ALLxsData)));
    end
    predfoo1 = mean(pred,2);
    predfoo1(predfoo1>0)=1;
    predfoo1(predfoo1<0)=-1;
    onemodel.trees=feattrees;
    onemodel.type='DF';
    result.featureTimeout=featLabel;
    result.predictfeatureTimeout=predfoo1;
    PBS.FeatModel=[];
    PBS.AllFeatLabel=predfoo1;
else
    PBS.FeatModel=[];
    result.featureTimeout=realFeatTime*0-1;
    result.predictfeatureTimeout=realFeatTime*0-1;
end

modelIDS=find(yDataOrig(:,result.Presolver(1))> result.PresolverTime(1) & yDataOrig(:,result.Presolver(2))>result.PresolverTime(2) & featLabel==-1);

modelIndex=1;

for jj = 1:result.NumSolver
    for kk = (jj+1):result.NumSolver
        fprintf('Build models for solver %d %d\n', jj, kk);
        trainfile=sprintf('%s.train_%d_%d', outfilename, jj, kk);
        modelfile=sprintf('%s.model_%d_%d', outfilename, jj, kk);
        testfile=sprintf('%s.test_%d_%d', outfilename, jj, kk);
        resultfile=sprintf('%s.result_%d_%d', outfilename, jj, kk);
        
        % first we need cost data and may be want somehow normalize the
        % cost data and we need label the data as -1 and 1
        % what is the cost
        traincostfoo=yDataSCORE(modelIDS,jj)-yDataSCORE(modelIDS,kk);
        trainlabelfoo=ones(length(traincostfoo), 1);
        trainlabelfoo(find(traincostfoo<=0))=-1;
        
        testcostfoo=ALLyDataSCORE(:,jj)-ALLyDataSCORE(:,kk);
        testlabelfoo=ones(length(testcostfoo), 1);
        testlabelfoo(find(testcostfoo<=0))=-1;
        
        traincostfoo=CostTransformation(traincostfoo, result.Cost);
        testcostfoo=CostTransformation(testcostfoo, result.Cost);
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
            fout = fopen(testfile, 'w');
            for uu=1:length(testlabelfoo)
                fprintf(fout, '%d cost:%f ', testlabelfoo(uu),  abs(testcostfoo(uu)));
                for ww=1:size(ALLxDataFinal,2)
                    fprintf(fout, '%d:%f ', ww, ALLxDataFinal(uu, ww));
                end
                fprintf(fout, '\n');
            end
            fclose(fout);
            
            mycmd=sprintf('!%ssvm_learn %s %s %s', svmpath, trainfile, modelfile);
            eval(mycmd);
            
            mycmd=sprintf('!%ssvm_classify %s %s %s', svmpath, testfile, modelfile, resultfile);
            eval(mycmd);
            predfoo2=csvread(resultfile);
            
            %% remove tmp file as we do that
            delete(trainfile);
            delete(testfile);
            delete(modelfile);
            delete(resultfile);
            onemodel.modelfile=modelfile;
            onemodel.type='SVM';
            onemodel.transformation=result.Transformation;
            PBS.Model{modelIndex}=[];
            PBS.AllPred{modelIndex}=predfoo2;
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
                X = xDataFinal(modelIDS(r),:);
                y = trainlabelfoo(r);
                cost = abs(traincostfoo(r));
                trees{num_tree} = classregtree(X ,y, 'method', 'classification', 'weights', cost, 'nvartosample', result.Nvartosample);
            end
            onemodel.trees=trees;
            onemodel.type='DF';
            onemodel.transformation=result.Transformation;
            
            pred_test=[];
            for num_tree=1:result.NumTree
                pred_test(:, num_tree) = str2num(char(eval(trees{num_tree}, ALLxDataFinal)));
            end
            predfoo2 = mean(pred_test,2);
            predfoo2(find(predfoo2>0))=1;
            predfoo2(find(predfoo2<0))=-1;
            
            PBS.Model{modelIndex}=[];
            PBS.AllPred{modelIndex}=predfoo2;
        end
        % no DT at this point
        %         if strcmp(result.Predictor, 'DT')
        %             t = classregtree(xDataFinal(modelIDS,:) ,trainlabelfoo, 'method', 'classification', 'weights', abs(traincostfoo));
        %             onemodel.tree=t;
        %             onemodel.type='DT';
        %             onemodel.transformation=result.Transformation;
        %             PBS.Model{modelIndex}=[];
        %         end
        result.PredictItem1{modelIndex}= jj;
        result.PredictItem2{modelIndex} = kk; 
        modelIndex=modelIndex+1;
    end
end
if result.NumSolver>1
    PBS.PredictItem1=result.PredictItem1;
    PBS.PredictItem2=result.PredictItem2;
else
    PBS.PredictItem1=[];
    PBS.PredictItem2=[];
end
