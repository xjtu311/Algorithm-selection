%% compute score for every canadidate solver based on some given reference
%% solvers. (The solvers which compete with SATzilla in SAT competition)
%% this is exactly time 2 score
function [yData, sScore] = time2score2(yDataAll, yReference, cutoff, useSCORE, serial, useSERIAL)
stdP = 1000;
spdP = 1000;
yData = yDataAll * 0;
[numI, numS]=size(yDataAll);
yDataAll = min(yDataAll, cutoff);
yReference = min(yReference, cutoff);
serM =3;
sScore =zeros(size(yDataAll));
if useSCORE ==1
    for ii= 1:numS;
        tmpYdata=[yDataAll(:,ii), yReference];
        foo = tmpYdata *0; ;
        foo(find(tmpYdata < cutoff)) =1;
        tmpScore = foo .* repmat(stdP ./(max(sum(foo,2), 0.0001)), 1, size(foo,2));
        spdf = cutoff ./ (1+tmpYdata);
        spdf(find(tmpYdata >=cutoff)) = 0;
        spds = spdP * spdf ./repmat(max(sum(spdf,2),0.0001), 1, size(foo,2));
        yData(:,ii) = tmpScore(:,1) + spds(:,1);
    end
    yData(find(yDataAll >=cutoff)) = 0;
    if useSERIAL ==1
        totalSerial = max(serial);
        for jj=1:numS
            tmpYdata=[yDataAll(:,jj), yReference];
            for ii = 1:totalSerial
                sID =find(serial == ii);
                sizeID =size(sID, 1);
                tmpsScore = zeros(sizeID,size(yReference,2)+1);
                if sizeID > 0
                    tmpyDataAll=tmpYdata(sID,:);
                    tmpsScore(find(tmpyDataAll<cutoff)) =1;
                    tmpSum = sum(tmpsScore,2);
                    tmpnum = sum(tmpsScore,1);
                    numGood = size(find(tmpnum>0),2); 
                    tmpSCORE = find(tmpSum > 0);
                    if size(tmpSCORE, 2) >0
                        if sizeID > 0
                            score = serM *stdP /(numGood+0.000001);
                        else
                            score = stdP /(0.000001+numGood);
                        end
                        tmpsScore = tmpsScore .* score;
                        sScore(min(find(tmpsScore(:,1)==1)),jj) = score;
                    end
                end
            end
        end
        sScore(find(yDataAll >=cutoff)) = 0;
    end
    
else
    yData = yDataAll;
end
