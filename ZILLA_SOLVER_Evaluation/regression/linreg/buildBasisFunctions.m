function [Phi, names, Aid, Bid] = buildBasisFunctions(X, namesX, square)

%=== Square parameters.
%[xParams, namesParams] = fh_squareData(xParams, namesParams);
%%[xParams, namesParams] = fh_squareData(xParams, namesParams);
%[xParams, namesParams, Aid, Bid] = squareData(xParams, namesParams);
fullX = X;
fullNames = namesX;

%numParams = size(xParams,2);

% pack; %stores stuff in contiguous blocks -> to counter memory problems ...
if square
    %=== Square everything (again).
    [Phi, names, Aid, Bid] = squareData(fullX, fullNames);
else
%     [Phi, names, Aid, Bid] = combineParamsAndFeats(fullX, fullNames, numParams);
   
   Phi = fullX;
   names = fullNames;
   Aid=1:length(fullNames);
   Bid=Aid*0;
end
%fprintf('Constructed %i basis functions.\n', size(Phi,2));