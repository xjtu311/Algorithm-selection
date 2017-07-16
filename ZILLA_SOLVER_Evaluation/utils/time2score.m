%% compute score for every canadidate solver based on some given reference
%% solvers. (The solvers which compete with SATzilla in SAT competition)

function yData = time2score(yDataAll, cutoff, scoreClass)
if strcmp(scoreClass, 'PAR10')    
    yData=yDataAll;
    yData(find(yData>=cutoff))=cutoff*10;
    yData=cutoff*10 - yData;
end

if strcmp(scoreClass, 'PAR1')    
    yData=yDataAll;
    yData(find(yData>=cutoff))=cutoff;
    yData=cutoff - yData;
end

if strcmp(scoreClass, 'LOG_PAR1')
    yData=log10(yDataAll);
   yData(find(yDataAll>=cutoff))=log10(cutoff);
   yData=log10(cutoff)-yData;

end

if strcmp(scoreClass, 'LOG_PAR10')
   yData=log10(yDataAll);
   yData(find(yDataAll>=cutoff))=log10(cutoff*10);
   yData=log10(cutoff*10)-yData;

end
