function [Y_mean, Y_var] = LRT(model, X)

[N,D] = size(X);
X = [X, ones(N,1)];

Y_mean = X * model.mu;
%=== Could get the whole covariance for the outputs if we wanted. 
%=== But this is an NxN matrix, big for lots of data points
%Y_var = diag(1/beta + X * Sigma_w * X');
Y_var = ones(size(Y_mean));