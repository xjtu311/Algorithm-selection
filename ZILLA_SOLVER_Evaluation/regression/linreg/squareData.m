function [sqX, sqNames, Aid, Bid] = squareData(X,names)

%code is similar to the one used in x2fx
[n m] = size(X);
first=repmat(1:m,m,1);
second=first';
t=first<=second;

sqX = zeros(n, m+m+m*(m-1)/2);
sqX = zeros(size(sqX,1), m+m+m*(m-1)/2);
Aid = zeros(1, m+m+m*(m-1)/2);
Bid = zeros(1, m+m+m*(m-1)/2);
sqX(:,1:m) = X;
Aid(:,1:m) = [1:m];
sqX(:,m+1:end) = X(:,first(t)) .* X(:,second(t));
Aid(:,m+1:end) = first(t)';
Bid(:,m+1:end) = second(t)';
sqNames = {};
if ~isempty(t)
    sqNames = strcat(names(first(t)),'.x.',names(second(t)));
end
sqNames = [names; sqNames];