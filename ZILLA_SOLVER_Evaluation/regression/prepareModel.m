function [model, nlml, dnlml] = prepareModel(model, invKmm)
assert(~any(isnan(model.y)));
assert(~any(isinf(model.y)));

switch model.type
    case {'rf', 'LR', 'ivm'}
    otherwise
        if isempty(strfind(model.type, 'GP'))
            error 'need to implement model for this modelType';
        end
end
    
if model.prepared
    nlml = 0;
    dnlml = 0;
    warning 'Model already prepared.'
%    return
end

switch model.type
    case 'LR'
        model.mu = LR(model.X, model.y);
        
    case 'ivm'
        model.ivm = learnIvmModel(model.X, model.y, model.cens);

    case 'regression-tree'
        %=== Normalize params to in between 0 and 1.
        params = (model.params - model.options.kernelParamLowerBound) / (model.options.kernelParamUpperBound-model.options.kernelParamLowerBound);
        split_min = 10;%ceil(50*params(1));  % 50 is the upper bound
        pruning = 'on';

        model.T = fh_treefit(model.X, model.y, 'splitmin', split_min, 'prune', pruning, 'catidx', model.cat);

    case 'rf'
        N = length(model.y);
        
        %=== Unpack params.
        ratioFeatures = model.params(1); 
        Splitmin = model.params(2);

        if any(model.cens)
            error 'The RF implementation does not handle censored data at this point.'
        else
            %=== Just do "normal" RF.
            for i=1:length(model.module)
                if model.options.orig_rf
                    r = ceil(N.*rand(N,1));
                else
                    r = 1:N; %new type of RF, where each tree has the same data. %ceil(N.*rand(perModelN,1));
                end
                Xsub = model.X(r,:);
                ysub = model.y(r);
                censsub = model.cens(r);

                seed = ceil(rand*10000000);
                model.module{i}.T = fh_simple_random_regtreefit_algoParams_instFeats_c(Xsub, ysub, censsub, Splitmin, ratioFeatures, model.cat, seed, model.options.logModel, model.catDomains);
            end
        end
        
end
    
if strfind(model.type, 'GP')
    global gprParams;
    gprParams = [];
    gprParams.combinedcat = model.cat;
    gprParams.combinedcont = model.cont;
    
    sigma_sqr_e = exp(2*model.params(end));
    model.var = exp(2*model.params(end-1));
    model.g = (model.var)/(model.var+sigma_sqr_e);
    
    if isfield(model.options, 'ppSize') && model.options.ppSize > 0
        %=== GP with subset of regressors or projected progress
        if any(model.cens)
            error('No censoring implemented yet for SRPP')
        end
        if nargout == 1
            if nargin > 1
                [model.Kmm, model.invKmm, model.saved_1, model.saved_2] = gprSRPPprepare(model.params, model.covfunc, model.X, model.pp_index, model.y, invKmm);
            else
                [model.Kmm, model.invKmm, model.saved_1, model.saved_2] = gprSRPPprepare(model.params, model.covfunc, model.X, model.pp_index, model.y); 
            end
        else
            error('prepareModel can only have 1 output if we use SRPP')
        end
    else
        %=== Normal GP
        if any(model.cens)
%        model.cens = zeros(length(model.y),1);
            model.useCensoring = 1;
            if nargout == 1
                [model.L_nonoise, model.alpha_nonoise, model.invK_times_invH_times_invK] = gprCensorPrepare(model.params, model.covfunc, model.X, model.y, model.cens);
            else
                error('prepareModel can only have 1 output if we use censoring')
            end
        else
            if nargout == 1
                [model.K, model.L, model.invL, model.alpha] = gprPrepare(model.params, model.covfunc, model.X, model.y);
            elseif nargout == 2
                [model.K, model.L, model.invL, model.alpha, nlml] = gprPrepare(model.params, model.covfunc, model.X, model.y);
            elseif nargout == 3
                [model.K, model.L, model.invL, model.alpha, nlml, dnlml, dKs] = gprPrepare(model.params, model.covfunc, model.X, model.y);
            else
                error('prepareModel needs to have at least 1 and not more than 3 outputs')
            end
        end
    end
end

model.prepared = 1;