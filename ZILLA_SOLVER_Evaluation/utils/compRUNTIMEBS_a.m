function [b,runtime] = compRUNTIMEBS(predv, realv, realcensor)
realv=min(realv, realcensor);
runtime = [];
numI=size(realv, 1);

[a, b]=max(predv, [], 2);
for i=1:numI
    runtime=[runtime; realv(i, b(i))];
end
