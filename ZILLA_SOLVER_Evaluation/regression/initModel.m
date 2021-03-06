function model = initModel(X, y, cens, cat, catDomains, isclean, options, names)
%=== Initialize and put data into model, but don't do any work yet.

if nargin < 7
    names = {};
    if nargin < 6
        options =  get_default_options();
    end
end

%=== Assert the inputs are fine.
assert(~any(isnan(y)));
assert(~any(isinf(y)));
assert(~any(isnan(cens)));
assert(~any(isinf(cens)));
assert(~any(any(isnan(X))));
assert(~any(any(isinf(X))));

%=== Clean the data: drop empty columns, normalize, PCA, quadratic basis functions, etc.
if options.logModel
    responseTransformation = 'log10';
else
    responseTransformation = 'identity';
end
linearSize = inf;
if isfield(options, 'linearSize')
    linearSize = options.linearSize;
end
maxModelSize = inf;
if isfield(options, 'maxModelSize')
    maxModelSize = options.maxModelSize;
end
doQuadratic = 0;
if isfield(options, 'doQuadratic')
    doQuadratic = options.doQuadratic;
end
numPcaComponents = 0;
if isfield(options, 'numPcaComponents')
    numPcaComponents = options.numPcaComponents;
end

perm = randperm(length(y));
trainIdx = perm(1:ceil(length(perm)/2));
validIdx = setdiff(1:length(perm), trainIdx);
xTrain = X(trainIdx,:);
yTrain = y(trainIdx);
censTrain = cens(trainIdx);
xValid = X(validIdx,:);
yValid = y(validIdx);
censValid = cens(validIdx);

if ~isclean
    [model.transformation, xOutput, yOutput, catOutput, catDomainsOutput] = buildCleanData(options, names, xTrain, yTrain, censTrain, xValid, yValid, censValid, cat, catDomains, numPcaComponents, linearSize, doQuadratic, maxModelSize, responseTransformation, 1);
    [xTmp, yTmp] = formatData(xTrain, yTrain, model.transformation);
    assert(max(max(abs(xTmp-xOutput))) < 0.00001);
    assert(max(abs(yTmp-yOutput)) < 0.00001);

    [X, y] = formatData(X, y, model.transformation);
    cat = catOutput;
    catDomains = catDomainsOutput;
end

%=== Remember inputs.
model.X = X;
model.y = y;
model.cens = cens;
model.cat = cat;
model.cont = setdiff(1:size(X,2), model.cat);
model.catDomains = {};
for i=1:length(catDomains)
    model.catDomains{i} = catDomains{i};
end
model.type = options.modelType;
model.options = options;
model.names = names;


%=== Model-specific initialization.

switch options.modelType
    case 'kriging'
        model.params = zeros(size(model.X,2)+1,1);
        model.params(end) = 0;
        %model.params(end) = options.kernelParamLowerBound + 0.99 * (options.kernelParamUpperBound-options.kernelParamLowerBound); %
        model.covfunc = {'covHybridSE_DIFFard'};

    case 'mixture'
        model.module = cell(1, options.nSub);
        opts = options;
        opts.modelType = opts.subModelType;
        tmpModel  = initModel(X, y, cens(1), model.cat, opts, model.names);
        model.params = tmpModel.params;
        %model.params = [tmpModel.params; 1.5]; % Extra parameter: fraction of data points to put in each tree.

    case 'rf'
        model.module = cell(1, options.nSub);
%        model.params = [];
        model.params = [options.init_split_ratio; options.init_Splitmin];
        %[options.regtree_p, options.split_ratio]; % pca (0.5*10)=5, p (10^(10*0.5-5) = 1), ratio feats. at each split (-1 = 1/3)
        %[(options.pca/10)*6-3, ((log10(options.regtree_p)+5)/10)*6-3, options.split_ratio*6-3]; % pca (0.5*10)=5, p (10^(10*0.5-5) = 1), ratio feats. at each split (-1 = 1/3)
        %        model.params = [-1.8,0,0]; % pca (0.5*10)=5, p (10^(10*0.5-5) = 1), ratio feats. at each split (0.5)
        

    case 'orig_bcm'
        use_normalization = 1;
        if use_normalization
            gpnet = gp(size(model.X,2), 'sqexp');
        else
            gpnet = gp(size(model.origX,2), 'sqexp');
        end
        net = bcm(gpnet);
        
        numClusters = ceil(size(model.origX,1)/options.subSize);%max(2, ceil(size(X,1)/options.subSize));
        
        kmeansopt = [1 1e-5 1e-4 0 0 0 0 0 0 0 0 0 0 30];
        r = randperm(size(model.X,1));
        if use_normalization
            [centres,opt,post] = kmeans(model.X(r(1:numClusters), :),model.X,kmeansopt);
        else
            [centres,opt,post] = kmeans(model.origX(r(1:numClusters), :),model.origX,kmeansopt);
        end
        [m,assignment] = max(post,[],2);
for i = 1:numClusters,
  fprintf('Clustered BCM: Module %i has %i data points\n', i, nnz(assignment==i));
end
        model.assignment = assignment;

        if use_normalization
            model.net = bcminit(net, model.X, model.origY, assignment);
        else
            model.net = bcminit(net, model.origX, model.origY, assignment);
        end
        model.params = '';
        
    case 'netlab'
        % Initialise the parameters.
        net = gp(size(model.X,2), 'sqexp'); %'ratquad'
%        prior.pr_mean = 0;
%        prior.pr_var = 1;
%        net = gpinit(net, x, t, prior);
        model.net = gpinit(net, model.X, y);
        model.params = gppak(model.net)';

    case 'LR'
        model.params = [];
        % [0.5,0.5]'; %delta (negative exponent divided by 10, i.e. 0=>10^0, 1=>10-10), max. subset size(multiplied by 50)
        
    case 'BLR'
        model.params = [0.5,0.5]'; %delta (negative exponent divided by 10, i.e. 0=>10^0, 1=>10-10), max. subset size(multiplied by 50)
        % for Bayesian linear regression
        [N,M] = size(model.X);
        model.beta = 1; % observation precision 1/sigma_obs^2.
        model.mean_0 = 0*ones(M+1,1); % Mx1 matrix, zero mean.
        model.Sigma_0 = 1e2 * eye(M+1,M+1); % MxM matrix, wide prior.

    case 'regression-tree'
        model.params = [0.5, 0.5]'; % order: split_min_factor (times 50), pruning (fraction of ntermnodes)

    case 'regtree'
        model.params = [-1.8];%[0]'; % order: split_min_factor (times 50)

    case 'randomregtree'
        model.params = [-3];  %default [-3; -1] [-1.8; 1];%[0]'; % order: split_min_factor (times 50)
        
    case 'GP-matern3'
        covfunc = {'covMatern3iso'};
        loghyper = zeros(2, 1);
    case 'GP-matern5'
        covfunc = {'covMatern5iso'};
        loghyper = zeros(2, 1);
    case 'GP-matern3noise'
        covfunc = {'covSum', {'covMatern3iso','covNoise'}};
        loghyper = zeros(3, 1);
    case 'GP-matern5noise'
        covfunc = {'covSum', {'covMatern5iso','covNoise'}};
        loghyper = zeros(3, 1);

    case 'GP-SEiso'
        covfunc = 'covSEiso';
        loghyper = zeros(2, 1);
    case 'GP-SEard'
        covfunc = {'covSEard'};
        loghyper = zeros(1 + size(model.X,2), 1);
    case 'GP-SEardnoise'
        covfunc = {'covSum', {'covSEard','covNoise'}};
        loghyper = zeros(2 + size(model.X,2), 1);
    case 'GP-SEisonoise'
        covfunc = {'covSum', {'covSEiso','covNoise'}};
        loghyper = zeros(3, 1);

    case 'GP-covRQardnoise'
        covfunc = {'covSum', {'covRQard','covNoise'}};
        loghyper = zeros(3 + size(model.X,2), 1);

    case 'GP-DIFFiso'
        covfunc = 'covDIFFiso';
        loghyper = zeros(2, 1);
    case 'GP-DIFFard'
        covfunc = 'covDIFFard';
        loghyper = zeros(1 + size(model.X,2), 1);

    case 'GP-hybridiso'
        covfunc = 'covHybridSE_DIFFiso';
        loghyper = zeros(2, 1);
    case 'GP-hybridard'
        covfunc = 'covHybridSE_DIFFard';
        loghyper = zeros(1 + size(model.X,2), 1);
    case 'GP-hybridisonoise'
        covfunc = {'covSum', {'covHybridSE_DIFFiso','covNoise'}};
        loghyper = zeros(3, 1);
    case {'GP-hybridardnoise', 'GPML'}
        covfunc = {'covSum', {'covHybridSE_DIFFard','covNoise'}};
        loghyper = zeros(2 + size(model.X,2), 1);
    case 'GP-hybridisoNoLen'
        covfunc = 'covHybridSE_DIFFiso_NoLen';
        loghyper = zeros(1, 1);
    case 'GP-hybridisoNoLennoise'
        covfunc = {'covSum', {'covHybridSE_DIFFiso_NoLen','covNoise'}};
        loghyper = zeros(2, 1);
    case 'GP-hybridisoOnlyLen'
        covfunc = 'covHybridSE_DIFFiso_OnlyLen';
        loghyper = zeros(1, 1);
    case 'GP-hybridisoOnlyLennoise'
        covfunc = {'covSum', {'covHybridSE_DIFFiso_OnlyLen','covNoise'}};
        loghyper = zeros(2, 1);

    case 'GP-covHybridSE_DIFFiso_algo'
        covfunc = {'covSum', {'covHybridSE_DIFFiso_algo','covNoise'}};
        loghyper = zeros(4, 1);

    case 'GP-AlarmMPE'
        covfunc = 'covAlarmMPE';
        loghyper = [log(1e-2)];
    case 'ogp'
        ogptrain(model.X,model.y);
    case 'dace'
        model.params = [];
    case 'dace-submodel'
        model.params = [];
        
    case 'dace-gpml'
        model.params = [];

    case 'ivm'
        model.params = [];
        
    otherwise
        error ('No such model type defined yet!');
end

%=== If we work with a subset of the data for optimization, set that index.
if isfield(options, 'trainSubSize') && options.trainSubSize > 0
    index = randperm(length(y));
    model.opt_index = index(1:min(length(y), options.trainSubSize));
end
%=== If we employ the projected process approximation, set the index for that.
if isfield(options, 'ppSize') && options.ppSize > 0
    index = randperm(length(y));
    model.pp_index = index(1:min(length(y), options.ppSize));
end

if strfind(options.modelType, 'GP')
    model.params = loghyper;

    if ischar(covfunc), covfunc = cellstr(covfunc); end % convert to cell if needed
    [n, D] = size(model.X);
    if eval(feval(covfunc{:})) ~= size(model.params, 1)
        error('Error: Number of parameters does not agree with covariance function')
    end
    model.covfunc = covfunc;

    if length(covfunc) == 2 & strcmp(covfunc(1), 'covSum') & strcmp(covfunc{2}(end), 'covNoise')
        model.noisy = 1;
        model.params(1:end-2) = 1; % start optimization with quite high length scale!
        model.params(end-1) = 0; % signal variance, average
        model.params(end) = -1; % start optimization with quite low noise!
    end
else
    model.options.ppSize = 0;
end

model.prepared = 0;

assert(~any(isnan(model.y)));
assert(~any(isinf(model.y)));