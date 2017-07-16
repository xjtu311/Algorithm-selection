%% let do something about option
if strcmp(result.Predictor, 'LR')
    options.modelType = 'LR';
    options.pca = 0;
    options.logModel = 0;
    options.doQuadratic = result.doQuadratic;
    options.linearSize = result.linearSize;
    options.maxModelSize = result.MaxModelSize;
end

if strcmp(result.Predictor, 'RF')
    options.modelType = 'rf';
    options.nSub = result.NumTree;
    options.logModel = 0;
    options.linearSize = 100000;
    options.maxModelSize = 100000;
    options.pca = 0; % do no PCA here, for feature importance
    options.init_split_ratio = result.init_split_ratio; % options.tuning_params(3);
    options.init_Splitmin = 10; % options.tuning_params(2);
    options.orig_rf = 0;
    options.strategyForMissing = 0;
    options.doQuadratic = 0;
end

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

for kk = 1:result.NumSolver
    metricfoo = time2score(yDataOrig(:,kk), result.Cutoff, result.MetricType);
    yDataMETRIC=[yDataMETRIC,metricfoo];
    metricfoo = time2score(ALLyData(:,kk), result.Cutoff, result.MetricType);
    ALLyDataMETRIC=[ALLyDataMETRIC,metricfoo];
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
    
    onemodel.trees=feattrees;
    onemodel.type='DF';
    
    PBS.FeatModel=[];
    pred=[];
    for num_tree=1:99
        pred(:, num_tree) = str2num(char(eval(feattrees{num_tree}, ALLxsData)));
    end
    predfoo1 = mean(pred,2);
    predfoo1(predfoo1>0)=1;
    predfoo1(predfoo1<0)=-1;
    result.featureTimeout=featLabel;
    result.predictfeatureTimeout=predfoo1;
    PBS.AllFeatLabel=predfoo1;
else
    PBS.FeatModel=[];
    result.featureTimeout=realFeatTime*0-1;
    result.predictfeatureTimeout=realFeatTime*0-1;
end

modelIDS=find(yDataOrig(:,result.Presolver(1))> result.PresolverTime(1) & yDataOrig(:,result.Presolver(2))>result.PresolverTime(2) & featLabel==-1);
namesX=[];
for fname=1:size(xDataOrig,2)
    sfoo=sprintf('%d', fname);
    namesX{fname}=sfoo;
end
namesX=namesX';
for ii = 1: result.NumSolver
    model = learnModel(xDataOrig(modelIDS,:), yDataSCORE(modelIDS,ii), yDataSCORE(modelIDS,ii)*0, [], [], 0, options, namesX);
    model.X=[];
    model.y=[];
    model.cens=[];
    model.names=[];
    tmpmodel=model;
    if strcmp(result.Predictor, 'RF')
        for hh=1:length(model.module)
            tmpmodel.module{1}.T.ysub=[];
            tmpmodel.module{1}.T.is_censored=[];
            tmpmodel.module{1}.T.nodeprob=[];
            tmpmodel.module{1}.T.risk=[];
            tmpmodel.module{1}.T.node=[];
            tmpmodel.module{1}.T.parent=[];           
            tmpmodel.module{1}.T.nodesize=[];
            tmpmodel.module{1}.T.npred=[];
            tmpmodel.module{1}.T.leaf_g=[];
            tmpmodel.module{1}.T.leaf_m=[];
            tmpmodel.module{1}.T.leaf_n=[];
            tmpmodel.module{1}.T.catcols=[];
            tmpmodel.module{1}.T.leaf_mean=[];
            tmpmodel.module{1}.T.leaf_var=[];
            tmpmodel.module{1}.T.emp_mean_at_leaf=[];
            tmpmodel.module{1}.T.nodeprob=[];
            tmpmodel.module{1}.T.minval=[];
            tmpmodel.module{1}.T.maxval=[];
            tmpmodel.module{1}.T.risk=[];  
        end
    end

    PBS.Model{ii}=[];
    result.PredictItem1(ii)= ii;
    result.PredictItem2(ii) = 0;
    linfoo= applyModel(model, ALLxData, 0, 0, 0);
    PBS.AllPred{ii}=linfoo;
    
end

PBS.PredictItem1=result.PredictItem1;
PBS.PredictItem2=result.PredictItem2;

