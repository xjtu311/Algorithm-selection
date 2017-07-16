%
% xx=[1, 2, 3; 2, 1, 4; 5, 2, 0; 3, 7, 9];
% yy= [1.3, 1.2, 2.3; 2.4, 4.1, 4.4; 5.2, 2.6, 3.0; 1.3, 3.7, 2.9];
%% return results is the configureration and expect score value
% method =1 use max, =0 use min
%% prediction is predicted score or runtime (method =1 using score otherwsie use runtime)
%% realVaule is the realruntime
%% referecne is the runtime for reference solvers
function [score,avgtime,solved] = solverSubsetSelectionSVMtest(tmpsolver, prediction, realValue, reference, serial, cutoff, scoreClass)
rand('state',sum(100*clock));
numSolver = size(realValue,2);
myindex=[];
tmpindex=[];
for e=1:numSolver
    for r=(e+1):numSolver
        myindex=[myindex, e*1000+r]; %#ok<AGROW>
    end
end
for g=1:length(tmpsolver)
    for h=(g+1):length(tmpsolver)
        tmpindex=[tmpindex, tmpsolver(g)*1000+tmpsolver(h)]; % the uniue index of solver combination/classification prediction
    end
end
[c, ia, ib] = intersect(tmpindex, myindex); % find which prediction we should use
tmpPrediction = prediction(:,ib);
tmpRealValue = realValue(:,tmpsolver); %find which solver we should use
score = zeros(size(realValue,1),1);

selected = selectSolverSVM(tmpPrediction, length(tmpsolver), 1); % pick solver using prediction
timeSpend =[];
for j=1:size(realValue,1)
    timeSpend=[timeSpend; tmpRealValue(j, selected(j))]; % find each the runtime of each instances
end

[scorefoo, avgtimefoo, solvedfoo]=time2scoretest(timeSpend, reference, cutoff, 1, serial, 1, scoreClass); % transfer to score
score=scorefoo(1);
avgtime=avgtimefoo(1);
solved=solvedfoo(1);