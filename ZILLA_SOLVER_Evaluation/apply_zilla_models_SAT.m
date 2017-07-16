function [crossTest, crossTrain]=apply_zilla_models_SAT(dataset, preid, seed, mymodel)
%% here mytype is training data or test data
if isdeployed
%     seed = str2double(seed);
%     preid = str2double(preid);
    mypath='/ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/';
else
    mypath='./';
end
PBS_option.FeatureCutoff=120;
PBS_option.MaxPreTime=0.05;

PBS_option.Presolver=[1,1];
PBS_option.PresolverTime=[0,0];
PBS_option.MiniRuntime=0.005;

PBS_option.NumCrossValidation=10;
PBS_option.Seed=seed;
PBS_option.Backupsolver=1;
PBS_option.NumTree=99;

PBS_option.init_split_ratio=[];
PBS_option.MetricType='PAR10';
PBS_option.Nvartosample=[];

%% linear regression
% PBS_option.Predictor='LR';
% PBS_option.MaxModelSize=20;
% PBS_option.doQuadratic=1;
% PBS_option.linearSize =30;
% PBS_option.ScoreType='LOG_PAR10';
% PBS_option.Cost='RAW';

%% RF
% PBS_option.init_split_ratio=[];
% PBS_option.NumTree=99;
% PBS_option.Predictor='RF';
% PBS_option.ScoreType='LOG_PAR10';
% PBS_option.Cost='RAW';

%% SVM
%  PBS_option.Predictor='SVM';
%  PBS_option.ScoreType='PAR10';
% PBS_option.Cost='RAW';

%% DF
% PBS_option.Predictor='DF';
% PBS_option.NumTree=33;
% PBS_option.Nvartosample=[];
% PBS_option.ScoreType='PAR10';
% PBS_option.Cost='RAW';
% PBS_option.Cost='UNIFORM';

%% new data for testing
rand('twister', seed);

myfile=sprintf('%sdata/%s',mypath,dataset);
inputfile=sprintf('%s.csv',myfile);
trainfile=sprintf('%s-train.csv',myfile);
testfile=sprintf('%s-test.csv',myfile);

% find a place for save all models
outdir=sprintf('%sSATMODEL', mypath);
dirname=sprintf('%s-%d-%d', dataset,seed,preid);

outdir=sprintf('%s/%s', outdir, dirname);

% this depends on different data sets
solved=3;
runtime=2;
if strcmp(dataset, 'SATHAND')
    numberSolver=15;
    % numberSolver=4;
    simfeatID=[1:2];
    featID=[1:115];
    feattimeID=[116:125];
    PBS_option.Cutoff=5000;
    mycheck=3; % check if instance is solved by presolver
end

if strcmp(dataset, 'SATRAND')
    numberSolver=9;
    % numberSolver=4;
    simfeatID=[1:2];
    featID=[1:115];
    feattimeID=[116:125];
    PBS_option.Cutoff=5000;
    mycheck=3; % check if instance is solved by presolver
end

if strcmp(dataset, 'SATINDU')
    numberSolver=18;
    % numberSolver=4;
    simfeatID=[1:2];
    featID=[1:115];
    feattimeID=[116:125];
    PBS_option.Cutoff=5000;
    mycheck=3; % check if instance is solved by presolver
end

if strcmp(dataset, 'CPAIOR')
    numberSolver=6;
    % numberSolver=4;
    simfeatID=[1:2];
    featID=[1:102];
    feattimeID=[103:106];
    PBS_option.Cutoff=3600;
    mycheck=1; % check if instance is solved by presolver
end
if strcmp(dataset, 'CRR')
    numberSolver=4;
    % numberSolver=4;
    simfeatID=[1:2];
    featID=[1:96];
    feattimeID=[97:100];
    PBS_option.Cutoff=3600;
    mycheck=1; % check if instance is solved by presolver
end

if strcmp(dataset, 'CPAIORNF')
    numberSolver=6;
    % numberSolver=4;
    simfeatID=[1:2];
    featID=[1:106];
    feattimeID=[107:111];
    PBS_option.Cutoff=3600;
    mycheck=1; % check if instance is solved by presolver
end
if strcmp(dataset, 'CRRNF')
    numberSolver=4;
    % numberSolver=4;
    simfeatID=[1:2];
    featID=[1:98];
    feattimeID=[99:103];
    PBS_option.Cutoff=3600;
    mycheck=1; % check if instance is solved by presolver
end
if strcmp(dataset, 'BMNF')
    numberSolver=4;
    % numberSolver=4;
    simfeatID=[1:2];
    featID=[1:108];
    feattimeID=[109:113];
    PBS_option.Cutoff=3600;
    mycheck=1; % check if instance is solved by presolver
end


alldata=csvread(inputfile, 1);
numOutputs=5;
allnames = textread(inputfile, '%s', 1, 'whitespace', '\n', 'bufsize', 10000);
allnames = strread(allnames{1},'%s','whitespace',',');
allnames = deblank(allnames);
namesY = allnames(2:numOutputs:numberSolver*5);
allnamesX = allnames(numberSolver*5+2:end);


solvedTable=zeros(size(alldata,1), numberSolver);
for i=1:numberSolver
    % what we call solved
    foo=solvedTable(:,i)+1;
    foo(alldata(:,i*5-5+solved)==0)=0;
    %     mean(foo)
    foo(alldata(:,i*5-5+runtime)>PBS_option.Cutoff-1)=0;
    alldata(alldata(:,i*5-5+solved)==0,i*5-5+runtime)=PBS_option.Cutoff+1;
    %     mean(foo)
    solvedTable(:,i)=foo;
end

% now, it is time to remove instances no one can solve and also
% preprocessor can solve
good=find(sum(solvedTable,2)>=1 & alldata(:,numberSolver*5+mycheck)>0);
bestsingle=max(sum(solvedTable));


fprintf('The precentage of instances solved by Oralce: %f\n',length(good)/size(solvedTable,1));
fprintf('The precentage of instances solved by best single solver: %f\n',bestsingle/length(good));

alldata=alldata(good,:);

permdata=alldata(randperm(length(good)),:);

ALLRT=permdata(:,[2:numOutputs:numberSolver*5]);
ALLSF=permdata(:,numberSolver*5+1+simfeatID);
ALLF=permdata(:,numberSolver*5+1+featID);
ALLfeattimefoo=permdata(:,numberSolver*5+1+feattimeID);
ALLfeattimefoo=max(ALLfeattimefoo,0);
ALLFT=sum(ALLfeattimefoo,2)*1;
%% new we can prepare data
startIdx=1;
testresult=[];
crossTest=[];
crossTest.AllFeatLabel=[];
trainresult=[];
crossTrain=[];
crossTrain.AllFeatLabel=[];

for l=1:PBS_option.NumCrossValidation
    % for l=1:1
    fprintf(strcat(['Cross-validation ', num2str(l), '/', num2str(PBS_option.NumCrossValidation), '**************************************\n']));
    endIdx = ceil(l*length(good)/PBS_option.NumCrossValidation);
    testIdx = startIdx:endIdx; % it is not the time to use test data yet
    trainIdx = setdiff(1:length(good), startIdx:endIdx);
    startIdx = endIdx+1;
    alltrain=permdata(trainIdx,:);
    RT=alltrain(:,[2:numOutputs:numberSolver*5]);
    SF=alltrain(:,numberSolver*5+1+simfeatID);
    F=alltrain(:,numberSolver*5+1+featID);
    feattimefoo=alltrain(:,numberSolver*5+1+feattimeID);
    feattimefoo=max(feattimefoo,0);
    FT=sum(feattimefoo,2)*1;
    
    %% with the same data, build models and save models for different settings
    
    %% linear regression
    if strcmp(mymodel,'LR')
        PBS_option.Predictor='LR';
        PBS_option.MaxModelSize=20;
        PBS_option.doQuadratic=1;
        PBS_option.linearSize =30;
        PBS_option.ScoreType='LOG_PAR10';
        PBS_option.Cost='RAW';
        outfilename=sprintf('%s-%s-%s-%d', PBS_option.Predictor,PBS_option.ScoreType,PBS_option.Cost,l);
        tmpresult=PBS_loader(RT, SF, F, FT, ALLRT, ALLSF, ALLF, ALLFT, PBS_option, outfilename,  outdir, preid);
        
    end
    %% RF
    if strcmp(mymodel,'RF')
        PBS_option.init_split_ratio=[];
        PBS_option.doQuadratic=0;
        PBS_option.linearSize =1000;
        PBS_option.MaxModelSize=1000000;
        PBS_option.NumTree=99;
        PBS_option.Predictor='RF';
        PBS_option.ScoreType='LOG_PAR10';
        PBS_option.Cost='RAW';
        outfilename=sprintf('%s-%s-%s-%d', PBS_option.Predictor,PBS_option.ScoreType,PBS_option.Cost,l);
        tmpresult=PBS_loader(RT, SF, F, FT, ALLRT, ALLSF, ALLF, ALLFT, PBS_option, outfilename,  outdir, preid);
        
    end
    %% SVM with cost
    if strcmp(mymodel,'SVM-RAW')
        PBS_option.Predictor='SVM';
        PBS_option.doQuadratic=0;
        PBS_option.linearSize =1000;
        PBS_option.MaxModelSize=1000000;
        PBS_option.ScoreType='PAR10';
        PBS_option.Cost='RAW';
        outfilename=sprintf('%s-%s-%s-%d', PBS_option.Predictor,PBS_option.ScoreType,PBS_option.Cost,l);
        tmpresult=PBS_loader(RT, SF, F, FT, ALLRT, ALLSF, ALLF, ALLFT, PBS_option, outfilename, outdir, preid);
        
    end
    %% SVM without cost
    if strcmp(mymodel,'SVM-UNI')
        PBS_option.Predictor='SVM';
        PBS_option.doQuadratic=0;
        PBS_option.linearSize =1000;
        PBS_option.MaxModelSize=1000000;
        PBS_option.ScoreType='PAR10';
        %     PBS_option.Cost='RAW';
        PBS_option.Cost='UNIFORM';
        outfilename=sprintf('%s-%s-%s-%d', PBS_option.Predictor,PBS_option.ScoreType,PBS_option.Cost,l);
        tmpresult=PBS_loader(RT, SF, F, FT, ALLRT, ALLSF, ALLF, ALLFT, PBS_option, outfilename,  outdir, preid);
        
    end
    %% DF with cost
    if strcmp(mymodel,'DF-RAW')
        PBS_option.Predictor='DF';
        PBS_option.NumTree=99;
        PBS_option.Nvartosample=[];
        PBS_option.doQuadratic=0;
        PBS_option.linearSize =1000;
        PBS_option.MaxModelSize=1000000;
        PBS_option.ScoreType='PAR10';
        PBS_option.Cost='RAW';
        % PBS_option.Cost='UNIFORM';
        outfilename=sprintf('%s-%s-%s-%d', PBS_option.Predictor,PBS_option.ScoreType,PBS_option.Cost,l);
        tmpresult=PBS_loader(RT, SF, F, FT, ALLRT, ALLSF, ALLF, ALLFT, PBS_option, outfilename,  outdir, preid);
       
    end
    %% DF without cost
    if strcmp(mymodel,'DF-UNI')
        PBS_option.Predictor='DF';
        PBS_option.NumTree=99;
        PBS_option.Nvartosample=[];
        PBS_option.doQuadratic=0;
        PBS_option.linearSize =1000;
        PBS_option.MaxModelSize=1000000;
        PBS_option.ScoreType='PAR10';
        % PBS_option.Cost='RAW';
        PBS_option.Cost='UNIFORM';
        outfilename=sprintf('%s-%s-%s-%d', PBS_option.Predictor,PBS_option.ScoreType,PBS_option.Cost,l);
        tmpresult=PBS_loader(RT, SF, F, FT, ALLRT, ALLSF, ALLF, ALLFT, PBS_option, outfilename,  outdir, preid);
    end
    %% for each iteration, we need remember the prediction
    testresult{l}.Presolver=tmpresult.Presolver;
    testresult{l}.PresolverTime=tmpresult.PresolverTime;
    testresult{l}.Predictor=tmpresult.Predictor;
    testresult{l}.Cutoff=tmpresult.Cutoff;
    testresult{l}.AllFeatLabel=tmpresult.AllFeatLabel(testIdx);
    for mm=1:length(tmpresult.AllPred);
        foo=tmpresult.AllPred{mm};
        testresult{l}.AllPred{mm}=foo(testIdx,:);
    end
    %% by the way, I also want know the orignal runtime data and feat runtime data
    testresult{l}.FT=ALLFT(testIdx,:);
    testresult{l}.RT=ALLRT(testIdx,:);
    %% Trainning performance
    trainresult{l}.Presolver=tmpresult.Presolver;
    trainresult{l}.PresolverTime=tmpresult.PresolverTime;
    trainresult{l}.Predictor=tmpresult.Predictor;
    trainresult{l}.Cutoff=tmpresult.Cutoff;
    trainresult{l}.AllFeatLabel=tmpresult.AllFeatLabel(trainIdx);
    for mm=1:length(tmpresult.AllPred);
        foo=tmpresult.AllPred{mm};
        trainresult{l}.AllPred{mm}=foo(trainIdx,:);
    end
    %% by the way, I also want know the orignal runtime data and feat runtime data
    trainresult{l}.FT=ALLFT(trainIdx,:);
    trainresult{l}.RT=ALLRT(trainIdx,:);
end
fprintf('Done loading the file!\n');
crossTest=testresult;
crossTrain=trainresult;
end
