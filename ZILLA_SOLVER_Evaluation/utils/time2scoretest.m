%% this code are used for compute the score for testing. it also compute
%% every reference score and so on
%% it use te same score function as before and treat the test data as a
%% real sat competition
%% the first one is the SATzilla runtime data
%% 10-7-2007

%% compute score for every canadidate solver based on some given reference
%% solvers. (The solvers which compete with SATzilla in SAT competition)

function [yData, avgTime, Solved] = time2scoretest(yDataAll, yReference, cutoff, scoreClass)
tmpYdata=[yDataAll, yReference,min(yReference,[],2)];
yDataAll = tmpYdata;
avgTime = mean(yDataAll,1);
if strcmp(scoreClass, 'PAR10')
    tmpYdata(find(tmpYdata>=cutoff))=cutoff*10;
    yData=cutoff*10 - tmpYdata;

end
if strcmp(scoreClass, 'PAR1')
    tmpYdata(find(tmpYdata>=cutoff))=cutoff*1;
    yData=cutoff*1 - tmpYdata;
end
if strcmp(scoreClass, 'LOG_PAR1')
    tmpYdata(find(tmpYdata>=cutoff))=cutoff*1;
    yData=log10(cutoff) - log10(tmpYdata);

end
if strcmp(scoreClass, 'LOG_PAR10')
    tmpYdata(find(tmpYdata>=cutoff))=cutoff*10;
    yData=log10(cutoff*10) - log10(tmpYdata);

end
solvedfoo=tmpYdata*0;
solvedfoo(find(tmpYdata<cutoff))=1;

Solved=sum(solvedfoo,1)/size(tmpYdata,1);
yData=sum(yData,1);
