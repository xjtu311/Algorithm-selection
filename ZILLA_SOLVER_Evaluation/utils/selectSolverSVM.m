%% This function pick the right solver given the classification results. if
%% the classification is perfect, the best one will have most vote.Even the
%% classification is not perfect, we can expect the best one have the most
%% votes


function pick=selectSolverSVM(pred, numsolver, mode)
%% add some random noise to the prediction to prevent ties
[afoo,bfoo]=size(pred);
pred=pred .* (ones(afoo, bfoo)+rand(afoo, bfoo)*0.0001);

pick=[];
pscore=pred*0;
nscore=pred*0;
pscore(find(pred>0))=1; %#ok<FNDSB>
nscore(find(pred<=0))=1;
switch mode
    case 1
        [instance, class]=size(pred);
        if class ~= numsolver*(numsolver-1)/2
            disp('Something Wrong with your prediction file\n');
        else
            vote=zeros(instance, numsolver);
            mygood=[];
            mybad=[];
            for ii=1:numsolver
                for jj=(ii+1):numsolver
                    mygood=[mygood, ii];
                    mybad=[mybad, jj];
                end
            end
            for tt=1:numsolver
               vote(:,tt)=sum(pscore(:, find(mygood==tt)),2)+sum(nscore(:,find(mybad==tt)),2);
            end
            %ok,now if there is a tie
            for mm=1:instance
                maxvote=max(vote(mm,:));
                maxsolver=find(vote(mm,:)==maxvote);
                if length(maxsolver)>1
%                     fprintf('There is a %d tie on instance %d!\n', length(maxsolver), mm);
                   foo=sort(maxsolver);
                   for aa=1:length(foo)
                       for bb= (aa+1):length(foo)
                           myfoo=find(mygood==foo(aa) & mybad==foo(bb));
                           myvalue=pred(mm,myfoo);
                           if myvalue>0
                               vote(mm, foo(aa))=vote(mm, foo(aa))+0.1+0.0001*abs(myvalue);
                           else
                              vote(mm, foo(bb))=vote(mm, foo(bb))+0.1+0.0001*abs(myvalue); 
                           end
                       end
                   end
                end
                % ok everything is good. Let's find the best one.
                %  pick=[pick; find(vote(mm,:)==max(vote(mm,:)))];
                bestone=find(vote(mm,:)==max(vote(mm,:)));
                pick=[pick; bestone(1)]; 
            end
        end
    otherwise
        fprintf('The mode %d is not supported yet!\n', mode);
end


end
