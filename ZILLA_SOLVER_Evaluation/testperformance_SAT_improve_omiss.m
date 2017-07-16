%% here is the code for findout presolver information and solver subset selection information

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



function alltestperformance=testperformance_SAT_improve_omiss(dataset, seed, mymodel, mode, omiss)
if isdeployed
    seed = str2double(seed);
    mode = str2double(mode);
    omiss = str2double(omiss);
end
mode=2; % always use mode 2

%% Now find out the right information
%% pay attenetion to missing data
allinfo=[];
for preid=1:144
    outputfilename=sprintf('./scripts-11-16-2011/SATzilla-SAT11-omiss/Train-%s-%s-%d-%d-%d.mat', dataset, mymodel, preid, seed, omiss);
    if exist(outputfilename)
        load(outputfilename);
        allinfo.preid{preid}.time=alltrainperformance.time{preid};
        allinfo.preid{preid}.score=alltrainperformance.score{preid};
        allinfo.preid{preid}.take=alltrainperformance.take{preid};
%        allinfo.preid{preid}=alltrainperformance;
    else
        allinfo.preid{preid}=[];
        fprintf('Missing preid %d data!\n', preid);
    end
end

%% now for each set of validation
prescore=zeros(144,10);
preidlist=[];
for cross=1:10
    for preid=1:144
        if ~isempty(allinfo.preid{preid})
            prescore(preid,cross)=allinfo.preid{preid}.score{cross}(1);
        end
    end
    bestfoo=find(prescore(:,cross)==max(prescore(:,cross)));
    preidlist=[preidlist,bestfoo(1)];
end

%% now find the best one (highest score)

% for a give data set
alltime=[];
allscore=[];
alltestperformance=[];
% for each cross validation
backupsolver=0;

for i=1:10
    [allmodeltest, allmodeldata]=apply_zilla_models_SAT_omiss(dataset, preidlist(i), seed, mymodel, omiss);
    modeldata=allmodeltest{i};
    numinst=length(modeldata.AllFeatLabel);
    onetime=zeros(numinst, 5);
    onesolve=onetime+1;
    pickfoo=[];
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
        alltestperformance.backupsolver{i}=backupsolver;
    else
      alltestperformance.backupsolver{i}=1; %only being used for INDU
      backupsolver=1;
    end
    onetime(modeldata.AllFeatLabel==1,1)=modeldata.RT(modeldata.AllFeatLabel==1,backupsolver);
    onesolve(modeldata.AllFeatLabel==1,[2:end])=0;
    alltestperformance.solvedbybackup{i}=find(modeldata.AllFeatLabel==1);

    % then deal with presolver 1
    solvedbypre1=find(modeldata.RT(:, modeldata.Presolver(1))<=modeldata.PresolverTime(1));
    onetime(:,2)=min(modeldata.RT(:, modeldata.Presolver(1)), modeldata.PresolverTime(1));
    onesolve(solvedbypre1,[3:end])=0;
    
    alltestperformance.Presolver{i}=modeldata.Presolver;
    alltestperformance.PresolverTime{i}=modeldata.PresolverTime;
    alltestperformance.solvedbypre1{i}=solvedbypre1;
    %pre2
    solvedbypre2=find(modeldata.RT(:, modeldata.Presolver(2))<=modeldata.PresolverTime(2));
    onetime(:,3)=min(modeldata.RT(:, modeldata.Presolver(2)), modeldata.PresolverTime(2));
    onesolve(solvedbypre2,[4:end])=0;

    alltestperformance.solvedbypre2{i}=solvedbypre2;
    
    % feature
    onetime(:,4)=modeldata.FT;
    alltestperformance.feattime{i}=modeldata.FT;
    % let's do model selection
    modelIDS=find(modeldata.RT(:,modeldata.Presolver(1))> modeldata.PresolverTime(1) & modeldata.RT(:,modeldata.Presolver(2))>modeldata.PresolverTime(2) & modeldata.AllFeatLabel==-1);
    alltestperformance.modelIDS{i}=modelIDS;
    validationTime=modeldata.RT+repmat(modeldata.FT, 1, size(modeldata.RT,2))+sum(modeldata.PresolverTime(1:2));
    allPred=[];
    for oo=1:length(modeldata.AllPred)
        allPred=[allPred, modeldata.AllPred{oo}];
    end
    if strcmp(mymodel, 'LR') || strcmp(mymodel, 'RF')
          PBStake=allinfo.preid{preidlist(i)}.take{i};
        [pickfoo, modeltime]=compRUNTIMEBS_a(allPred(:,PBStake),  modeldata.RT(:,PBStake), modeldata.Cutoff);
    else
        PBStake=allinfo.preid{preidlist(i)}.take{i};
        [pickfoo, modeltime]=compRUNTIMESVM_a(PBStake, allPred,  modeldata.RT);
    end
    fprintf('finish crossvalidation %d\n', i);
    
    onetime(:,5)=modeltime;
    testperformance=sum(onetime .* onesolve, 2);
    
    [PBSscore, PBStime,  PBSsolved] =time2scoretest(testperformance, modeldata.RT, modeldata.Cutoff, 'PAR10');
    
    alltestperformance.picksolver{i}=pickfoo;
    alltestperformance.time{i}=testperformance;
    alltestperformance.score{i}=PBSscore;
    alltestperformance.take{i}=PBStake;
end
if mode==2
outputfilename=sprintf('Test1-%s-%s-%d-%d', dataset, mymodel, seed, omiss);
end
save(outputfilename, 'alltestperformance');


