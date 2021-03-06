function model = optimizeModel(model, varargin)
%=== Optimize the model, e.g. by optimizing kernel parameters.
options = model.options;
if ~isfield(options, 'opt') || ~options.opt
    return
end

if strfind(model.type, 'dace')
    return; % optimization done in prepare.
end

if model.options.trainSubSize > 0
    %=== Optimize the model parameters on a subset of the data.
    opts = options;
    opts.trainSubSize = 0;
    opts.ppSize = 0; % Don't use PP for the optimization.
%    subModel = initModel(model.origX(model.opt_index,:), model.origY(model.opt_index), model.cens(model.opt_index), model.cat, opts, model.origNames);
    subModel = subsetModel(model, model.opt_index);
    subModel.options = opts;
    subModel = optimizeModel(subModel);
    
    %=== Then prepare a model with those hyperparameters.
    model.params = subModel.params;
    return
end

if isempty(model.params)
    return
end

numHypParams = length(model.params);
%==== Different objectives for hyperparameter optimization.
switch options.hyp_opt_obj
    case 'mll' % marginal likelihood
        if strcmp(model.type, 'BCM') || strcmp(model.type, 'mixture')
            func_handle = @(params) summedNmllModel(params, model);
        else
            func_handle = @(params) nmllModel(params, model);
%            func_handle = @(params) nmllModelFixed_SF_and_noise(params, model);
        end
    case 'cv-ll' % cross-validated log-likelihood of unseen data.
%            if strcmp(model.type, 'BCM') || strcmp(model.type, 'mixture')
%                func_handle = @(params) summedCrossValLL(params, model);
%            else
            func_handle = @(params) crossValLL(params, model);
%            end
    case 'cv-rmse' % cross-validated RMSE of unseen data.
%            if strcmp(model.type, 'BCM') || strcmp(model.type, 'mixture')
%                func_handle = @(params) summedCrossValLL(params, model);
%            else
            func_handle = @(params) crossValRMSE(params, model);
%            end
    case {'marg', 'marg-rmse', 'marg-cc'}
        func_handle = @(params) validateModelWithParams(params, model, options.hyp_opt_obj, varargin{:});
    otherwise
        error strcat(['Unknown objective for hyperparameter optimization', options.hyp_opt_obj]);
end

%==== Different approaches for hyperparameter optimization.
switch options.hyp_opt_algorithm
    case 'minFunc'
        %=== minFunc setup.
        minFuncOptions.Method = 'lbfgs';
        minFuncOptions.MaxIter = options.hyp_opt_steps; % 50;
        minFuncOptions.numDiff = 0; % for now since I don't have a derivative
        if strcmp(model.type, 'kriging')
            minFuncOptions.numDiff = 1;
        end
        minFuncOptions.useComplex = 0;
        minFuncOptions.Display = 'on';
        minFuncOptions.DerivativeCheck = 'off';

        start_params = model.params;
        
        [params, obj] = minFuncBC(func_handle, start_params, options.kernelParamLowerBound*ones(numHypParams,1), (options.kernelParamUpperBound-1e-10)*ones(numHypParams,1), minFuncOptions);

    case 'cma-es'
        %=== CMA-ES setup.
        cmaes_opts = cmaes('defaults');
        cmaes_opts.Plotting = 'off';
        cmaes_opts.LBounds = options.kernelParamLowerBound;
        cmaes_opts.UBounds = options.kernelParamUpperBound;
        % Don't need sigma if we have bounds for each param       sigma = 1/3 * (options.kernelParamUpperBound-options.kernelParamLowerBound);

        %TODO: should #function evals by options.hyp_opt_steps
        
        for i=1:options.hyp_opt_numTries
            if i==1
                start_params = model.params;
            else
                start_params = randn(numHypParams,1);
            end

            switch options.hyp_opt_obj
                case 'mll'
                    [paramChoices(:,i), objs(i)] = cmaes('nmllModel', start_params, [], cmaes_opts, model);
                case 'cv-ll'
                    [paramChoices(:,i), objs(i)] = cmaes('crossValLL', start_params, [], cmaes_opts, model);
            end
        end
        objs
        [minObj, idx] = min(objs);
        params = paramChoices(:, idx);

    case 'direct'
        %=== DIRECT setup.
        opts.maxevals = options.hyp_opt_steps;
        opts.maxits = 1e10; % deactivated
        opts.showits = 1;
        bounds(1:numHypParams,1) = options.kernelParamLowerBound;
        bounds(1:numHypParams,2) = options.kernelParamUpperBound;
        Problem.f = func_handle;
        [obj, params, hist] = Direct(Problem, bounds, opts);

        %=== Output to log.
        global log_fid
        if ~isempty(log_fid)
            fprintf(log_fid, '%d\t%d\t%f\n', hist', params', obj');
        end
    otherwise
        error(strcat(['Unknown optimizer for kernel hyperparameters: ', options.hyp_opt_algorithm]));
end
model.params = params;

%=== HELPER FUNCTION TO GET THE GRADIENT FOR SCG
function res = second_output_as_row(f,x)
[tmp, res] = f(x);
if size(res,1) > 1 && size(res,2) > 1
    error 'Gradient has to be a vector';
end
res = res(:)';