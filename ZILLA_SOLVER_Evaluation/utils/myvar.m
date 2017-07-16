function a=myvar(x)
n = length(x);
a= 1/(n-1.0) * (sum(x.^2) - (1.0/n)*sum(x).^2);
