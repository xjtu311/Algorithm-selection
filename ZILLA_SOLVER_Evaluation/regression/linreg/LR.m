function mu= trainModel(trainX, trainY)

%Train regression model and return coefficients

error(nargchk(2,2,nargin));
delta=1e-2;
% delta=10;
TrainM = [trainX ones(size(trainX,1),1)];

% deal with rank deficiency, approximates using pseudo-inverse
A = TrainM' * TrainM + delta * eye(size(TrainM,2));
% condNumber=klbcond(A);

mu = inv(A) * (TrainM' * trainY);



