function [sqX, sqNames] = combineParamsAndFeats(X,names,numParms)

newX = [];
newNames = {};
%X=X(1:10,:);
[N,D] = size(X);
numFeats = D-numParms;
for i=1:numParms
    newX = [newX, repmat(X(:,i),[1,numFeats]) .* X(:,numParms+1:end)];
    if nargout>1
        tmpNames = strcat(strcat(names{i}, '.x.'),names(numParms+1:end));
        newNames = [newNames; tmpNames];
    end
end

sqX = [X,newX];
if nargout>1
    sqNames = [names; newNames];
end