function count = checkSubsetCell(cacheOfMeanAndVar,cacheOfMeanAndVar2)
count = 0;
for i=1:length(cacheOfMeanAndVar)
    if ~isempty(cacheOfMeanAndVar{i})
        tmpCell = cacheOfMeanAndVar{i};
        for j=1:length(cacheOfMeanAndVar)
            if ~isempty(tmpCell{j})
                tmpCell2 = tmpCell{j};
                for k=1:length(tmpCell2)
                    assert( length(cacheOfMeanAndVar{i}{j}{k}) == length(cacheOfMeanAndVar2{i}{j}{k}) );
                    for l=1:length(cacheOfMeanAndVar{i}{j}{k})
                        assertVectorEq(cacheOfMeanAndVar{i}{j}{k}{l}, cacheOfMeanAndVar2{i}{j}{k}{l}, 1e-4);
                    end
                    if ~isempty(tmpCell2{k})
                        count = count+1;
                    end
                end
            end
        end
    end
end
               
                        
                        
                
        
    
    