%% Portfolio based algorithm selection for HAL
% input: runtime
%        features
%        featuretime
%        options for PBS
%        output file
% output: PBS information in output file
%         some information in the results
%% example:  ob.Predictor='LR'; ob.Cutoff=3600; ob.ScoreType='PAR-10';
%           aa=PBS_builder(rand(100,3), rand(100,3), rand(100, 1), ob, 'ss');
%% At this stage, this is a pure test function. It has no value except for
%% debugging purpose.
function result=PBS_loader(RT, SF, F, FT, ALLRT, ALLSF, ALLF, ALLFT, PBS_option, filename, outputdir, preid)

% some constence such as svm path
tic; % time information
tStart = tic;

% svmpath='./svm/';
if isempty(PBS_option.Presolver)
    result.Presolver=[1,1];
    result.PresolverTime=[0.0, 0.0];
else
    result.Presolver=PBS_option.Presolver;
    result.PresolverTime=PBS_option.PresolverTime;
end

if isempty(PBS_option.Predictor)
    result.Predcitor='DF';
    result.NumTree = 99;
    result.ScoreType = 'PAR10';
    result.Cost = 'SQRT';
    result.Metric = 'PAR10';
else
    result.Predictor=PBS_option.Predictor;
    result.NumTree = ceil(PBS_option.NumTree/2)*2-1; %Make sure it is an odd number
    result.MaxModelSize = PBS_option.MaxModelSize;
    result.ScoreType = PBS_option.ScoreType;
    result.Cost = PBS_option.Cost;
    result.MetricType = PBS_option.MetricType;
    result.doQuadratic = PBS_option.doQuadratic;
    result.linearSize = PBS_option.linearSize;
end

if isempty(PBS_option.init_split_ratio)
    result.init_split_ratio=5.0/6;
else
    result.init_split_ratio = PBS_option.init_split_ratio;
end

if isempty(PBS_option.MiniRuntime)
    result.MiniRuntime=0.005;
else
    result.MiniRuntime = PBS_option.MiniRuntime;
end

if isempty(PBS_option.NumCrossValidation)
    result.NumCrossValidation = 10;
else
    result.NumCrossValidation = PBS_option.NumCrossValidation;
end

if isempty(PBS_option.Seed)
    result.Seed=1234;
else
    result.Seed=PBS_option.Seed;
end

result.NumSolver=size(RT,2);
result.NumInstance=size(RT,1);
result.NumFeature=size(F,2);
result.NumSimpleFeature=size(SF,2);
result.MaxPreTime=PBS_option.MaxPreTime;

result.Cutoff=PBS_option.Cutoff;
result.FeatureCutoff=PBS_option.FeatureCutoff;

if result.FeatureCutoff <=0
    if isempty(PBS_option.Backupsolver)
        tmpavg=mean(min(max(RT,result.MiniRuntime), result.Cutoff));
        result.Backupsolver=find(tmpavg==min(tmpavg));
        result.Backupsolver = result.Backupsolver(1); % in case, we have ties.
    else
        result.Backupsolver=PBS_option.Backupsolver;
    end
else
    if isempty(PBS_option.Backupsolver)
        tmp=min(max(RT,result.MiniRuntime), result.Cutoff);
        tmp1=tmp(sum(FT,2) <0);
        if isempty(tmp1)
            tmpavg=mean(min(max(RT,result.MiniRuntime), result.Cutoff));
            result.Backupsolver=find(tmpavg==min(tmpavg));
            result.Backupsolver = result.Backupsolver(1); % in case, we have ties.
        else
            tmpavg=mean(tmp1);
            result.Backupsolver=find(tmpavg==min(tmpavg));
            result.Backupsolver = result.Backupsolver(1); % in case, we have ties.
        end
    else
        result.Backupsolver=PBS_option.Backupsolver;
    end
end


% let do something about option
% train model to make preidction (it should be easy and we only want the best model. we do not do cross validation)
featnamesX=[];
for fname=1:size(SF,2)
    sfoo=sprintf('%d', fname);
    featnamesX{fname}=sfoo;
end
featnamesX=featnamesX';
result.FeatNvartosample = ceil(log2(size(SF,2))+1); % default for classification forest

%% get data already have permutation
% rand('state', result.Seed);
% perm = randperm(result.NumInstance);
yData = RT;
xData = F;
xTime = FT;
xsData=SF;
yData= min(max(yData, result.MiniRuntime), result.Cutoff);


ALLyData = ALLRT;
ALLxData = ALLF;
ALLxsData=ALLSF;
ALLxTime = ALLFT;
ALLyData= min(max(ALLyData, result.MiniRuntime), result.Cutoff);

% feature labels (not really going to use it)
realFeatTime=sum(xTime,2);
realFeatTime(realFeatTime<0)=result.FeatureCutoff*10;
featLabel=realFeatTime*0-1;
featLabel(realFeatTime>=result.FeatureCutoff)=1;

% backup copy of all the  data
yDataOrig = yData;
xDataOrig = xData;
xTimeOrig = xTime;

% record runing information
runtimeTable=zeros(size(xData,1), 5);
solvedTable=runtimeTable+1;
% let's decide what kind of presolver and pretime
if preid>0
    minavg=min(mean(yDataOrig,1));
    maxpretime= max(minavg * result.MaxPreTime, 0.5);
    maxpresolved=[];
    for i=[1:result.NumSolver]
        maxpresolved=[maxpresolved, length(find(yDataOrig(:,i) <= maxpretime))];
    end
    pretimes=[0, maxpretime/4.0, maxpretime/2.0, maxpretime];
    [a,b]=sort(maxpresolved*-1);
    presolvers=b(1:min(length(b),3));
    
    pretables=[];
    for a=1:length(presolvers)
        for b=1:length(presolvers)
            for c=1:4
                for d=1:4
                    pretables=[pretables;[presolvers(a),presolvers(b),pretimes(c),pretimes(d)]];
                end
            end
        end
    end
    result.Presolver=[pretables(preid, 1),pretables(preid, 2)];
    result.PresolverTime=[pretables(preid, 3), pretables(preid, 4)];
end

%% I should be able to take something out for the combination of presolvers
% if I decided not to run that, we should be able tell which one is the
% equvient
todolist=ones(size(pretables,1),1);

for jjj=1:length(todolist)
    haveit=0;
    comptodolist=find(todolist==0);
    for iii=1:length(comptodolist)
        mytodo=pretables(jjj,:);
        thistodo=pretables(comptodolist(iii),:);
        % pre1=pre2', pre2=pre1', pret1=pre2t', pret2=pre1t'
        if mytodo(1)==thistodo(2) && mytodo(2)==thistodo(1) ...
                && mytodo(3)==thistodo(4) && mytodo(4)==thistodo(3)
            haveit=comptodolist(iii);
        end
        % 0 only have once
        if haveit==0
            if  mytodo(3)== 0 &&mytodo(4)== 0 && mytodo(3)==thistodo(4) &&mytodo(4)==thistodo(3)
                haveit=comptodolist(iii);
            end
        end
        % pre1=pre2
        if haveit==0
             if  mytodo(1)== mytodo(2)
                 realcutoff=max(mytodo(3),mytodo(4));
                 if mytodo(1)==thistodo(1) && thistodo(3)==realcutoff && thistodo(4)==0
                     haveit=comptodolist(iii);
                 end
                  if mytodo(1)==thistodo(2) && thistodo(4)==realcutoff && thistodo(3)==0
                     haveit=comptodolist(iii);
                 end
             end
        end
        if haveit==0
            if  mytodo(3)== 0
                
                if mytodo(2)==thistodo(1) && thistodo(3)==mytodo(4) && thistodo(4)==0
                    haveit=comptodolist(iii);
                end
                if mytodo(2)==thistodo(2) && thistodo(4)==mytodo(4) && thistodo(3)==0
                    haveit=comptodolist(iii);
                end
            end
        end 
        if haveit==0
            if  mytodo(4)== 0
                
                if mytodo(1)==thistodo(1) && thistodo(3)==mytodo(3) && thistodo(4)==0
                    haveit=comptodolist(iii);
                end
                if mytodo(1)==thistodo(2) && thistodo(4)==mytodo(3) && thistodo(3)==0
                    haveit=comptodolist(iii);
                end
            end
        end 
    end
    todolist(jjj)=haveit;
end


if todolist(preid)>0
    %% we need modify output dir
     outputdir=regexprep(outputdir, '\d*$', int2str(todolist(preid)));
    
end
global svmpath;
outfilename=sprintf('%s/%s', outputdir,filename);
svmpath='/ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/svm/';

%% we need rememeber all the result stuff now
load(outfilename);
tmpresult=PBS;
tmpresult.Presolver=result.Presolver;
tmpresult.PresolverTime=result.PresolverTime;
tmpresult.Cutoff=result.Cutoff;
result=tmpresult;

% [foo,invperm]=sort(perm);
% trainSolverTime=runtimeTable(invperm, 5);
% solvedTablefoo=solvedTable(invperm,5);
% trainSolverTime(find(solvedTablefoo==0))=-512;


% result.PBSAvgRuntime=PBStime(1);
% result.OracleAvgRuntime=PBStime(end);
% result.SolverAvgRuntime=PBStime(2:end-1);
% result.FeatAvgRuntime=mean(FT);
% result.PBSScore=PBSscore(1);
% result.OracleScore=PBSscore(end);
% result.SolverScore=PBSscore(2:end-1);
% result.PBSSolvedPercentage=PBSsolved(1);
% result.OracleSolvedPercentage=PBSsolved(end);
% result.SolverSolvedPercentage=PBSsolved(2:end-1);
% result.SolverPicked=PBStake;
% result.Elapsed = toc(tStart);
% result.TrainSolverTime=trainSolverTime;
% PBS.SolverPicked=PBStake;

% in order to make tar work. I need to change the working dir. Then change
% it back
% currpath=pwd;
% chdir(outputdir);
% tarfilename=sprintf('tar-%s.tgz', filename);
% targetfile=sprintf('%s.*', filename);
% mycmd=sprintf('!tar -zcvf %s %s', tarfilename, targetfile);
% eval(mycmd);
% unix(['rm -f ' targetfile]);
% chdir(currpath);

end