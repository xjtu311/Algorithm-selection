%
% xx=[1, 2, 3; 2, 1, 4; 5, 2, 0; 3, 7, 9];
% yy= [1.3, 1.2, 2.3; 2.4, 4.1, 4.4; 5.2, 2.6, 3.0; 1.3, 3.7, 2.9];
%% return results is the configureration and expect score value
% method =1 use max, =0 use min
%% prediction is predicted score or runtime (method =1 using score otherwsie use runtime)
%% realVaule is the realruntime
%% referecne is the runtime for reference solvers
function output = solverSubsetSelectionSVM(prediction, realValue, cutoff, scoreClass)
rand('state',sum(100*clock));
numSolver = size(realValue,2);
restart = 20;
randness =0.95;
stopped=400;
results =[];
configure =[];
stdP = 1000;
spdP = 1000;
pick =zeros(1, numSolver);
ccheck =0;
myindex=[];

for e=1:numSolver
    for r=(e+1):numSolver
        myindex=[myindex, e*1000+r]; %#ok<AGROW>
    end
end

for i =1:restart
    % rand pick some solvers and put into pick
    while sum(pick) ==0
        pick = rand(1, numSolver);
        pick(find(pick>0.5))=1;
        pick(find(pick <= 0.5)) = 0;
    end
    bbest=pick*0+1;
    abest = 0;
    ccheck =1;
    acounter=0;
    cont = 1;
    while cont==1
        %% the prediction is more complicated in the case of using one-to-one classfication
        tmpsolver=find(pick==1); % picked solvers
        tmpindex=[];
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
        
        [scorefoo,avgtime, solved]=time2scoretest(timeSpend, [], cutoff, scoreClass); % transfer to score
        score=scorefoo(1);

        if abest <= score % if it is better, take it
            if (abest == score) && (sum(bbest) <= sum(pick))
            else
                abest = score;
                tmpPick =pick;
                bbest = pick;
                ccheck = 1;
            end
        end

        if(rand(1)>randness) % if it is not bad, still take it with randness but not change the abest and bbest
            tmpPick=pick;
        else
            pick =tmpPick;
            [foo, cid] = max(rand(1, numSolver));
            pick(cid)=abs(pick(cid)-1);
            while(sum(pick)==0) %at least take one
                [foo, cid] = max(rand(1, numSolver));
                pick(cid)=abs(pick(cid)-1);
            end
        end
        ccheck =ccheck+1;
        if ccheck > stopped
            cont =0;
        end
        acounter = acounter +1;
        %        fprintf('The best score so far is %f \n', tmpResult);
    end
    %         fprintf('The best score so far is %f with counter %d\n', abest, acounter);
    results=[results; abest];
    configure = [configure; bbest];
end


[a,b] = sort(results,1, 'descend');

output.configure = configure(b,:);
output.value = results(b);



