%% this is the code for do local search based sovler subset selection to obtain the best oracle performace
% input is a set of performance such as runtime, par, and the size of the
% set
% output is the set itself and the oracle performance with that set

function output=picksubsetsolvers(data, maxsize)
rand('state',sum(100*clock));
numSolver = size(data,2);
restart = 40;
randness =0.9;
stopped=10000;
results =[];
configure =[];
ccheck =0;
for i =1:restart
    pick =zeros(1, numSolver);
    init=randperm(numSolver);
    init=init(1:maxsize);
    pick(init)=1;
    
    abest = 1.0e99;
    bbest=zeros(1, numSolver);
    
    ccheck =1;
    acounter=0;
    cont =1;
    
    while cont==1
        score=sum(min(data(:,pick>0), [],2));
        % if get better performance take it
        if abest > score
            abest = score;
            bbest = pick;
            ccheck = 1;
            
            % now modify pick to get a new pick
            turnoff=randperm(maxsize);
            onitem=find(pick==1);
            pick(onitem(turnoff(1)))=0;
            turnon=randperm(numSolver-maxsize);
            offitem=find(pick==0);
            pick(offitem(turnon(1)))=1;
            
            
        else
            if(rand(1)>randness)
                % in some cases, we want keep this one even it is not the best
                % now modify pick to get a new pick
                turnoff=randperm(maxsize);
                onitem=find(pick==1);
                pick(onitem(turnoff(1)))=0;
                turnon=randperm(numSolver-maxsize);
                offitem=find(pick==0);
                pick(offitem(turnon(1)))=1;
                
                
            else
                ccheck =ccheck+1;
                if ccheck > stopped
                    cont =0;
                end
                acounter = acounter +1;
                %        fprintf('The best score so far is %f \n', tmpResult);
            end
            %         fprintf('The best score so far is %f with counter %d\n', abest, acounter);
            
        end
        
        
        
    end
    results=[results; abest];
    configure = [configure; bbest];
end
[a,b] = sort(results);

output.configure = configure(b,:);
output.value = results(b);
end
