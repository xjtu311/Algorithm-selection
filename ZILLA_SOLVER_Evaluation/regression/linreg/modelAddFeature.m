% function [model A Ainv XtY] = modelAddFeature(X, new_feature, y, A, Ainv, delta, XtY)
% Train the linear model (X'X + I \delta)^{-1}X'y, given the matrix product
% and inverse from the previous iteration
function [model,A,Ainv,XtY]  = modelAddFeature(X, new_feature, y, A, Ainv, delta, XtY)
	%X = [X new_feature];
	%temp = X(:,1:end-1)' * X(:,end);
	%A = [A temp; temp' X(:,end)'*X(:,end)];  
	%persistent tempX aBar aDot u1 v1 u2 v2;
	tempX = X' * new_feature;
	A = [A tempX; tempX' new_feature'*new_feature];  
	A(end,end) = A(end,end)+delta;

	% update the inverse stored in Ainv and train model
	p = size(A,1)-1;

	% augment the matrix (is there a better way?)
	Ainv(:,p+1) = zeros;
	Ainv(p+1,:) = zeros;
	Ainv(p+1,p+1) = 1;

	aBar = A(p+1,1:p); % row
	aDot = A(p+1,p+1);
	
	% first rank one update
	u1 = [zeros(p,1); 1];
	v1 = [aBar'; 0];
	%temp =  v1'*Ainv;
	%Ainv = Ainv - (Ainv(:,end)*temp)/(1 + temp(:,end));
	Ainv = Ainv - ((Ainv*u1)*(v1'*Ainv))/(1 + v1'*Ainv*u1);
	
	% second rank one update
	u2 = [aBar'; aDot-1];
	v2 = u1;
	%Ainv = Ainv - ((Ainv*u2)*(v2'*Ainv))/(1 + v2'*Ainv*u2);
	Ainv = Ainv - ((Ainv*u2)*(v2'*Ainv))/(1 + v2'*Ainv*u2);
	
	% finally, build the model
	XtY = [XtY; new_feature' * y];
	model = Ainv * XtY;
