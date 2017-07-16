%% compute score first
yDataSCORE=[];
testyDataSCORE=[];
validationtime=[];

for kk = 1:(nSOLVER)
    [scorefoo, foo2] = time2score(yDataORIG(:,kk), yDataORIG(:, reference), cutoff, 1, ...
        zeros(size(xDataORIG,1),1), useSERIAL, 5);
    validationtime=[validationtime, testyDataORIG(:,kk)+testxDataORIG(:,end)+(preTime1+preTime2)];
    yDataSCORE=[yDataSCORE,scorefoo];
    [scorefoo, foo2] = time2score(testyDataORIG(:,kk), testyDataORIG(:, reference), cutoff, 1, ...
        zeros(size(testxDataORIG,1),1), useSERIAL, 5);
    testyDataSCORE=[testyDataSCORE,scorefoo];
end

transformation =normalData(xDataORIG(:, goodfeatures));
for jj = 1:nSOLVER
    xDataFinal=formatnormalData(xDataORIG(:, goodfeatures),transformation);
    testxDataFinal=formatnormalData(testxDataORIG(:, goodfeatures),transformation);
end
%% all predictins, we will compute RMSE later
result.TRAINPRED=[];
result.TESTPRED=[];
modelIDSAll=[];

%% have cross validation
startIdx = 1;
for l=1:ck
    N=length(perm);
    fprintf(strcat(['Cross-validation ', num2str(l), '/', num2str(ck), '**************************************\n']));
    endIdx = ceil(l*N/ck);
    testIdx = startIdx:endIdx;
    trainIdx = setdiff(1:N, startIdx:endIdx);
    startIdx = endIdx+1;

    %% first do presolving. Now we need to have cross validation

    modelIDS=trainIdx(find(yDataORIG(trainIdx, preSolver1)>preTime1 & yDataORIG(trainIdx, preSolver2)>preTime2));
    modelIDSAll=[modelIDSAll, testIdx(find(testyDataORIG(testIdx, preSolver1)>preTime1 & testyDataORIG(testIdx, preSolver2)>preTime2))];
    %% for test data we need update runtimeTable and solvedTable

    %pre1

    solvedbypre1=testIdx(find(testyDataORIG(testIdx, preSolver1)<=preTime1));
    runtimeTable(testIdx,1)=min(testyDataORIG(testIdx, preSolver1), preTime1);
    solvedTable(solvedbypre1,[2:end])=0;

    %pre2
    solvedbypre2=testIdx(find(testyDataORIG(testIdx, preSolver2)<=preTime2));
    runtimeTable(testIdx,2)=min(testyDataORIG(testIdx, preSolver2), preTime2);
    solvedTable(solvedbypre2,[3:end])=0;

    % feature
    runtimeTable(testIdx,3)=testxData(testIdx,end);
    tmpTrainPred=[];
    tmpTestPred=[];
    %totally. 2) treat them as 0 in the normalized data.


    for jj = 1:nSOLVER
        for kk = (jj+1):nSOLVER
            trainfile=sprintf('%strain%f', svmpath, rand());
            testfile=sprintf('%stest%f', svmpath, rand());
            modelfile=sprintf('%smodel%f', svmpath, rand());
            resultfile=sprintf('%sresult%f', svmpath, rand());
            % first we need cost data and may be want somehow normalize the
            % cost data and we need label the data as -1 and 1
            % what is the cost
            testcostfoo=testyDataSCORE(testIdx,jj)-testyDataSCORE(testIdx,kk);
            testlabelfoo=ones(length(testcostfoo), 1);
            testlabelfoo(find(testcostfoo<=0))=-1;
            traincostfoo=yDataSCORE(modelIDS,jj)-yDataSCORE(modelIDS,kk);
            trainlabelfoo=ones(length(traincostfoo), 1);
            trainlabelfoo(find(traincostfoo<=0))=-1;

            if useWeight == 0
                testcostfoo=testcostfoo*0+1;
                traincostfoo=traincostfoo*0+1;
            end
            if useWeight == 1
                testcostfoo=abs(testcostfoo);
                traincostfoo=abs(traincostfoo);
            end
            if useWeight ==2
                traincostfoo=sqrt(abs(traincostfoo)); %% put less weight
                testcostfoo=sqrt(abs(testcostfoo));
            end

            fprintf(strcat(['Cross-validation ', num2str(l), '/', num2str(ck), '**************************************\n']));

            if useSVM==1
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
                    for ww=1:size(testxDataFinal,2)
                        fprintf(fout, '%d:%f ', ww, testxDataFinal(testIdx(uu), ww));
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

            if (useSVM ==2)
                num_trees = 99;
                numInst = length(modelIDS);
                nvartosample = ceil(log2(size(xDataFinal,2))+1); % default for classification forest
                for num_tree=1:num_trees
                    r = ceil(numInst.*rand(numInst,1));
                    X = xDataFinal(modelIDS(r),:);
                    y = trainlabelfoo(r);
                    cost = abs(traincostfoo(r));
                    trees{num_tree} = classregtree(X ,y, 'method', 'classification', 'weights', cost, 'nvartosample', nvartosample);
                end
                pred=[];
                pred_test=[];
                for num_tree=1:num_trees
                    pred(:, num_tree) = str2num(char(eval(trees{num_tree}, xDataFinal(modelIDS,:))));
                    pred_test(:, num_tree) = str2num(char(eval(trees{num_tree}, testxDataFinal(testIdx,:))));
                end

                predfoo1 = mean(pred,2);
                predfoo2 = mean(pred_test,2);
                predfoo1(find(predfoo1>0))=1;
                predfoo1(find(predfoo1<0))=-1;
                if (length(find(predfoo1==0)) > 0)
                    error('Tie. Use an odd number of trees to avoid this.')
                end

                predfoo2(find(predfoo2>0))=1;
                predfoo2(find(predfoo2<0))=-1;
                if (length(find(predfoo2==0)) > 0)
                    error('Tie. Use an odd number of trees to avoid this.')
                end
            end

            if useSVM==0
                t = classregtree(xDataFinal(modelIDS,:) ,trainlabelfoo, 'method', 'classification', 'weights', abs(traincostfoo));
                predclass1=eval(t, xDataFinal(modelIDS,:));
                predclass2=eval(t, testxDataFinal(testIdx,:));
                predfoo1=str2num(char(predclass1));
                predfoo2=str2num(char(predclass2));
            end

            tmpTrainPred=[tmpTrainPred, predfoo1];
            tmpTestPred=[tmpTestPred, predfoo2];
        end
    end
    result.TRAINPRED=[result.TRAINPRED; tmpTrainPred];
    result.TESTPRED=[result.TESTPRED; tmpTestPred];

end
%% let's see the crossvalidation result
indfoo=1;
for ii = 1: nSOLVER
    for jj = (ii+1): nSOLVER
        loss=compLossSVM(testyDataSCORE(:,ii), testyDataSCORE(:,jj), result.TESTPRED(:,indfoo));
        fprintf('Realtive loss for classification between %s and %s is: %f\n', namesY{ii}, namesY{jj}, loss);
        indfoo=indfoo+1;
    end
end

%% With all the predictions, we want know which solver we need to pick
%% based on validation performance
% the input will be predictions of which one is better n*(n-1)/2 of them
% output the solver should be use
if it>1
    result.ordered = solverSubsetSelectionSVM(result.TESTPRED(modelIDSAll,:), validationtime(modelIDSAll,alltake),  testyDataORIG(modelIDSAll,:), ...
        zeros(size(modelIDSAll, 1),1), cutoff, 5);

    takefoo=find(result.ordered.configure(1,:)==1);
else
    takefoo = [1];
end
result.TESTMODEL=compRUNTIMESVM(takefoo, result.TESTPRED,  testyDataORIG);

runtimeTable(:,4)=result.TESTMODEL;
result.TEST=sum(runtimeTable .* solvedTable, 2);

[score.TESTSCORE1, avgTime.TEST,  solved.TEST] =time2scoretest(result.TEST,   testyDataORIG(:, reference), cutoff, 1,zeros(size(result.TEST, 1),1),1,5);
for k=2:size(score.TESTSCORE1, 2)
    fprintf('Avg time for only using solver %s is %f and score is %f, and solved precentage is %f \n', ...
        namesY{k-1}, avgTime.TEST(k), score.TESTSCORE1(k), solved.TEST(k));
end
%% we also want output oracle result
oracleruntime=min(testyDataORIG,[], 2);

fprintf('Avg time for only using Oralce is %f and solved precentage is %f \n', ...
    mean(oracleruntime), length(find(oracleruntime<cutoff))/length(oracleruntime));


fprintf('Avg time for testing is %f, %f and score is %f, and solved predentage is %f*****\n',...
    avgTime.TEST(1), mean(result.TEST), score.TESTSCORE1(1), solved.TEST(1));

fprintf('Taked solvers are: ');
for i=1:length(takefoo)
    fprintf('%d, ', takefoo(i));
end
fprintf('\n');


fprintf('RESULT: %s, %d, %d, %f, %d, %f, %f, %f, %f\n', dataset, nSOLVER, preSolver1, preTime1, preSolver2, preTime2,  avgTime.TEST(1), mean(result.TEST),cutoff);

%% now I want par score for each instance instead of each run. What can we
%% do? We need to know every solvers and for oracle as well
allruntime=[result.TEST, testyDataORIG];
allmedian=[];
allpar=[];
for hh=1:size(allruntime,2)
    myfoo=allruntime(:,hh);
    allmedian=[allmedian, median(myfoo)'];
    myfoo(find(myfoo>=cutoff))=cutoff*10;
    allpar=[allpar, mean(myfoo)'];
end
