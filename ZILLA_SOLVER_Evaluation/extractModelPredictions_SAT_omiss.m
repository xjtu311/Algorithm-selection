%% this is the code for collect all the data for presolver selection
% The input will be a data set
%                   a model
% The output will be the result including presolver selection and
% solver subset selection

% dataset='RAND';
% dataset='HAND';
% dataset='INDU';
% dataset='CPAIOR';
% dataset='CRR';
% mymodel='SVM-RAW'
% mymodel='LR'
% mymodel='RF'
% mymodel='SVM-UNI'
% mymodel='SVM-RAW'
% mymodel='DF-RAW'
% mymodel='DF-UNI'



function alltrainperformance=extractModelPredictions_SAT_omiss(dataset, preid, seed, mymodel, omiss)
if isdeployed
    seed = str2double(seed);
    preid = str2double(preid);
    omiss = str2double(omiss);
end

% for a give data set
alltime=[];
allscore=[];
alltrainperformance=[];
% for each cross validation

% for each pre-id find the best validation performance
fprintf('Preid is %d, seed is %d\n', preid, seed);
[allmodeltest, allmodeldata]=apply_zilla_models_SAT_omiss(dataset, preid, seed, mymodel, omiss);
for i=1:10
    modeldata=allmodeldata{i};
    numinst=length(modeldata.AllFeatLabel);
    onetime=zeros(numinst, 5);
    onesolve=onetime+1;
    % well, we have to deal with backup solver first
    % who has the best performance on training data only appliable for INDU
    if strcmp(dataset, 'SATINDU')
        if omiss <11
          backupsolver=10;
        end
        if omiss >11
          backupsolver=11;
        end
        if omiss ==11
          backupsolver=11; % solver 12 is the backup solver
        end
    else
      backupsolver=1; %only being used for INDU
    end
    onetime(modeldata.AllFeatLabel==1,1)=modeldata.RT(modeldata.AllFeatLabel==1,backupsolver);
    onesolve(modeldata.AllFeatLabel==1,[2:end])=0;
    % then deal with presolver 1
    solvedbypre1=find(modeldata.RT(:, modeldata.Presolver(1))<=modeldata.PresolverTime(1));
    onetime(:,2)=min(modeldata.RT(:, modeldata.Presolver(1)), modeldata.PresolverTime(1));
    onesolve(solvedbypre1,[3:end])=0;
    
    %pre2
    solvedbypre2=find(modeldata.RT(:, modeldata.Presolver(2))<=modeldata.PresolverTime(2));
    onetime(:,3)=min(modeldata.RT(:, modeldata.Presolver(2)), modeldata.PresolverTime(2));
    onesolve(solvedbypre2,[4:end])=0;
    
    % feature
    onetime(:,4)=modeldata.FT;
    
    % let's do model selection
    modelIDS=find(modeldata.RT(:,modeldata.Presolver(1))> modeldata.PresolverTime(1) & modeldata.RT(:,modeldata.Presolver(2))>modeldata.PresolverTime(2) & modeldata.AllFeatLabel==-1);
    validationTime=modeldata.RT+repmat(modeldata.FT, 1, size(modeldata.RT,2))+sum(modeldata.PresolverTime(1:2));
    allPred=[];
    for oo=1:length(modeldata.AllPred)
        allPred=[allPred, modeldata.AllPred{oo}];
    end
    if strcmp(mymodel, 'LR') || strcmp(mymodel, 'RF')
        modeldata.ordered = solverSubsetSelectionLSR(allPred(modelIDS,:), validationTime(modelIDS,:),  modeldata.Cutoff, 'PAR10');
        PBStake=find(modeldata.ordered.configure(1,:)==1);
        
        modeltime=compRUNTIMEBS(allPred(:,PBStake),  modeldata.RT(:,PBStake), modeldata.Cutoff);
    else
        modeldata.ordered = solverSubsetSelectionSVM(allPred(modelIDS,:), validationTime(modelIDS,:),  modeldata.Cutoff, 'PAR10');
        PBStake=find(modeldata.ordered.configure(1,:)==1);
        
        modeltime=compRUNTIMESVM(PBStake, allPred,  modeldata.RT);
    end
    fprintf('finish crossvalidation %d\n', i);
    
    onetime(:,5)=modeltime;
    trainperformance=sum(onetime .* onesolve, 2);
    
    [PBSscore, PBStime,  PBSsolved] =time2scoretest(trainperformance, modeldata.RT, modeldata.Cutoff, 'PAR10');
    
    alltime{i}=trainperformance;
    allscore{i}=PBSscore;
    allpick{i}=PBStake;
end
alltrainperformance.time{preid}=alltime;
alltrainperformance.score{preid}=allscore;
alltrainperformance.take{preid}=allpick;
outputfilename=sprintf('Train-%s-%s-%d-%d-%d', dataset, mymodel, preid, seed, omiss);
save(outputfilename, 'alltrainperformance');

