%
% xx=[1, 2, 3; 2, 1, 4; 5, 2, 0; 3, 7, 9];
% yy= [1.3, 1.2, 2.3; 2.4, 4.1, 4.4; 5.2, 2.6, 3.0; 1.3, 3.7, 2.9];
%% return results is the configureration and expect score value
% method =1 use max, =0 use min
%% prediction is predicted score or runtime (method =1 using score otherwsie use runtime)
%% realVaule is the realruntime
%% referecne is the runtime for reference solvers
function output = solverSubsetSelectionLSR(prediction, realValue, cutoff, scoreClass)
rand('state',sum(100*clock));
numSolver = size(prediction,2);
restart = 10;
randness =0.95;
stopped=100;
results =[];
configure =[];
stdP = 1000;
spdP = 1000;
pick =zeros(1, numSolver);
ccheck =0;
for i =1:restart
    while sum(pick) ==0
        pick = rand(1, numSolver);
        pick(find(pick>0.5))=1;
        pick(find(pick <= 0.5)) = 0;
    end
    abest = 0;
    bbest=zeros(1, numSolver);
    tmpPick=zeros(1, numSolver);
    ccheck =1;
    acounter=0;
    cont =1;
    while cont==1
        tmpPrediction = prediction(:,find(pick==1));
        tmpRealValue = realValue(:,find(pick==1));
        score = zeros(size(prediction,1),1);

        [a, b]=max(tmpPrediction, [], 2);
        timeSpend =[];
        for j=1:size(prediction,1)
            timeSpend=[timeSpend; tmpRealValue(j, b(j))];
        end
        timeSpend=timeSpend;
        [scorefoo,avgtime, solved]=time2scoretest(timeSpend, [], cutoff, scoreClass);
        score=scorefoo(1);
        if abest <= score
            if (abest == score) & (sum(bbest) <= sum(pick))
            else
                abest = score;
                tmpPick =pick;
                bbest = pick;
                ccheck = 1;
            end
        end

        if(rand(1)>randness)
            tmpPick=pick;
        else

            pick =tmpPick;
            [foo, cid] = max(rand(1, numSolver));
            pick(cid)=abs(pick(cid)-1);
            while(sum(pick)==0)
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



