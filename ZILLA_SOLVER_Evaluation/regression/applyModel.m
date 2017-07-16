function [complete_yPred, complete_yPredVar, samples_of_min] = applyModel(model, X, isclean, observationNoise, joint)
if ~isclean
    if ~isfield(model, 'transformation')
        error('ERROR. The model was learnt on clean data, it doesnt have a transformation for unclean data.');
    end
    X = formatData(X, [], model.transformation);
end

if isfield(model, 'special_for_trivial') && (model.special_for_trivial==1)
    nontrivial_insts = [];
    for i=1:size(X,1)
        if ~any(X(i,:)==-512)
            nontrivial_insts = [nontrivial_insts, i];
        end
    end
    origSizeOfTest = size(X,1);
    X = X(nontrivial_insts,:);
end
yPred = -ones(size(X,1),1);
yPredVar = -ones(size(X,1),1);

if nargin < 4
    if nargin < 3
        observationNoise = 0;
    end
    joint = 0;
end

switch model.type
    case 'ivm'
        [yPred, yPredVar] = ivmActualPosteriorMeanVar(model.ivm, X);
    case 'BCM'
        querySize = model.options.subSize;
        verbosity = 0;
        [yPred, yPredVar] = fh_bcmFwd(model,X,querySize,verbosity);

    case 'mixture'
        yPredModels = zeros(size(X,1), length(model.module));
        for m=1:length(model.module)
            yPredModels(:,m) = applyModel(model.module{m}, X);
        end
        yPred = mean(yPredModels,2);
%        yPred = trafo(mean(inv_trafo(yPredModels, model.options.transformation),2), model.options.transformation);

        yPredVar = std(yPredModels,0,2).^2;
        smallIdx = find(yPredVar < 1e-10);
        yPredVar(smallIdx) = 1;
        yPredVar(smallIdx) = min(yPredVar);
        yPred = yPred + model.meanToAdd;

%        return
    case 'kriging'
        %p = [0.024894131360485201, 0.12500139809322541, 1.0002318133068488, 5.1658529179130808];
        %g = 0.28110580391601903;
%         p = [0.066666675319801680];
%         g = 0.99999363323199497;
%         model.params = [log(sqrt(1./p))'; 1]

%         g = 0.50728481462785591;
%         p =  [0.020412830067804148         0.12500130071390392        4.3157848670263128        10.207871395803409];
%         model.params = [log(sqrt(1./p))'; 1];

        g = model.g;
        [Kss, Kstar] = covHybridSE_DIFFard([model.params(1:end-1); 0], model.X, X); % 0 in log space => length scale 1
        Kstar = Kstar .^ 2;
        Kstar = Kstar .* g;

%         % Can normally just reuse the old R; but since we manually change
%         % the parameters for debugging, need to recompute.
%         R = covHybridSE_DIFFard([model.params(1:end-1); 0], model.X); % 0 in log space => length scale 1
%         R = R .^ 2;
%         R = R .* g; % diag should be one before this! last param is g
%         R = R - diag(diag(R)) + diag(diag(ones(length(R))));
%         
%         
%         
%         
%         %%%%%%%%%%%%%%tmp
%         
%         
%         [U,D] = eig(R);
%         invD = inv(D);
%         eigvalfit = diag(D); % eigs
%         diaginvfit = 1./eigvalfit; % inverse eigenvalues
% 
%         % From SKO: 
%         % compute beta hat as given by the mle				
%         %   beta  = (F' Rinv F)^(-1) (F' Rinv y)		
%         %         = (F' U Dinv U' F)^(-1) (F' U Dinv U' y)		
%         % 	 = ( F' U Dinv (F' U)' )^(-1)  ( F' U Dinv (y' U)' )
%         F = ones(length(U), 1);
%         tmp = F' * U;
% %        beta = inv( tmp * invD * tmp' ) * ( tmp * invD * (model.y' * U)' );
%         beta = inv( tmp .* diaginvfit' * tmp' ) * ( tmp .* diaginvfit' * (model.y' * U)' );
% 
%         % From SKO: 
%         % var         = 1/n ( (y-mu)' invR (y-mu) )			
%         %             = 1/n [U'(y-F beta)]' D^(-1) [U'(y-F beta)] 	
%         %             = 1/n [(y-F beta)' U] D^(-1) [(y-F beta)' U]'
%         % FH: special case of F=ones
% 
%         tmp = (model.y-beta)' * U;
%         var = 1/length(model.y) *  tmp * invD * tmp';
%         
%         V = R.*var;
%         
%         
%         %%%%%%%%%%%%%
        
        
%     Raugm = [0, ones(1, size(model.X,1)); ones(size(model.X,1),1 ), V];
%     invRaugm = inv(Raugm);
% 
%     for i=1:size(X,1)
%         yPred(i) = beta + Kstar(:,i)' * inv(V) * (model.y-beta);
%         
%         vec1 = [1; Kstar(:,i)];
%         yPredVar(i) = var*g - vec1' * invRaugm * vec1;        
%     end    


        Raugm = [0, ones(1, size(model.X,1)); ones(size(model.X,1),1 ), model.R];
        [Uaugm,Daugm] = eig(Raugm); % returns different order than the SKO implementation
        inv_Daugm = inv(Daugm);
        
        for i=1:size(X,1)
            vec1 = [1; Kstar(:,i)];
            yvec = [0; model.y];

            % From SKO: compute temp1 = (vec1' U) Diag(eigval)^{-1}.   need
            % for both MSE and BLUP
            temp1 = (vec1' * Uaugm) .* diag(inv_Daugm)';

            % From SKO: compute BLUP
            % multiply blup  = (vec1' U) Diag(eigval)^{-1}* (yvec' U)'
            % multiply       = temp1                      *  temp2     
            yPred(i,1) = temp1 * (yvec' * Uaugm)';

            % From SKO: compute MSE 
            % multiply s    = (vec1' U) Diag(eigval)^{-1}*  (vec1' U)'
            %               = temp1                      *     temp2  
            s = temp1 * (vec1' * Uaugm)';
            yPredVar(i,1) = model.var * (1.0 -  s); % includes observation noise
        end
        
        if ~observationNoise
            yPredVar = yPredVar - (model.var - g*model.var);
        end
        
    case 'rf'
% Do this in the calling function.
%         if model.pca > 0
%             %=== Apply PCA to instance features.
%             X_inst = X(:,model.numThetaLeft+1:end);
%             X_inst = X_inst(:,model.sub);
%             X_inst = X_inst - repmat(model.means, [N,1]);
%             X_inst = X_inst./repmat(model.stds, [N,1]);
%             X_inst = X_inst*model.pcVec;
%             X = [X(:,1:model.numThetaLeft), X_inst];
%         end
        
        samples_of_min = inf * ones(1, length(model.module));
        numchunks = ceil(size(X,1)*length(model.module) / 1000000);
        lenchunks = ceil(size(X,1)/numchunks);
        lower = 1;
        for i=1:numchunks
            upper = min([lower+lenchunks, size(X,1)]);
            subX = X(lower:upper,:);

            %=== Regular prediction for subX.
            treemeans = zeros(size(subX,1), length(model.module));
%           treemeans2= zeros(size(subX,1), length(model.module));
            for m=1:length(model.module)
                treemeans(:,m) = fh_simple_one_treeval(model.module{m}.T, subX, model.options.strategyForMissing);
                assert(~any(isnan(treemeans(:,m))));
                assert(~any(isinf(treemeans(:,m))));
%                 Tree = model.module{m}.T;
%                 cell_of_leaves = fh_treeval_thetas_pis(Tree, X, []);
%                 for i=1:length(cell_of_leaves)
%                     leaf_mean = Tree.emp_mean_at_leaf(cell_of_leaves{i}{4});
%                     treemeans2(cell_of_leaves{i}{1},m) = leaf_mean;
%                 end
%                 assertVectorEq(treemeans2(:,m), treemeans(:,m));                
            end
%                 assertVectorEq(treemeans, treemeans2);

%             if model.options.logModel
%                 subyPred = mean(log10(treemeans),2);
%                 subyPredVar = var(log10(treemeans),0,2);
%                 subsamples_of_min = min(log10(treemeans),[],1);
%             else
                subyPred = mean(treemeans,2);
                subyPredVar = var(treemeans,0,2);
                subsamples_of_min = min(treemeans,[],1);
%             end
            
            %=== Regular prediction for subX.
            yPred(lower:upper,1) = subyPred;
            yPredVar(lower:upper,1) = subyPredVar;
            samples_of_min = min([samples_of_min;subsamples_of_min],[],1);
            lower = upper+1;
            if lower > size(X,1)
                break;
            end
        end
        
        assert(~any(isnan(yPredVar)));
        assert(~any(isinf(yPredVar)));
        
    case {'dace'}
%         for i=1:size(X,1)
%             [yPred(i), dy0, yPredVar(i)] = predictor(X(i,:), model.subModel.dacemodel);
%         end
%         yPred = yPred(:);
%         yPredVar = yPredVar(:);
%         samples_of_min = min(model.merged_y);
        [yPred, yPredVar, samples_of_min] = applyModel(model.subModel, X, observationNoise, joint); 
        if any(yPredVar<0)
            minYPredVar = min(yPredVar)
            warning('yPredVar < 0');
            yPredVar = max(yPredVar, 1e-10);
        end
            

    case 'dace-submodel'
%         for i=1:size(X,1)
%             [yPred(i), dy0, yPredVar(i)] = predictor(X(i,:), model.dacemodel);
%         end
        if size(X,1) == 1
            [yPred, dy0, yPredVar] = predictor(X, model.dacemodel);
        else
            [yPred, yPredVar] = predictor(X, model.dacemodel);
        end
        toc
        yPred = yPred(:);
        yPredVar = yPredVar(:);
        samples_of_min = min(model.y);
        
    case 'dace-gpml'
        [yPred, yPredVar] = gprFwd(model.subModel.X, model.subModel.L, model.subModel.invL, model.subModel.alpha, model.subModel.covfunc, model.subModel.params, X, 0, 0);        
        samples_of_min = min(model.merged_y);
        
    case 'gp_eta'
        
        % setup...
        [nx, na] = size(model.X);
        uno = ones(nx, 1);
        delta2inv = inv(model.delta2);

        % first, compute r
%        r = computeCorrFast(model.X, model.theta, 0, X);
        r = ones(nx,1);
        for j = 1:nx
           r(j) = exp(-sum((1/model.theta) * (model.X(j,:)-X).^2));
        end

        % now predict the response
        yPred = model.mu + r' * model.invR * (model.y - model.uno * model.mu);

        % and the variance
        s = model.sig2 * (1 - r' * model.invR * r + (1 - model.uno' * model.invR *r)^2 / (model.uno' * model.invR* model.uno + inv(model.delta2)));
        yPredVar = s;
        
        yPred = yPred + model.meanToAdd;

        %return;

    case 'orig_bcm'
%        [yPred, yPredVar] = bcmfwd(model.net, X, 500); % not using the normalization
        [yPred, yPredVar] = bcmfwd(model.net, X, 500); % using the normalization

    case 'netlab'
        if joint
            K11 = gpcovarp(model.net, X, X);
            K12 = gpcovarp(model.net, X, model.X);
            % Prediction
            yPred = K12*model.net.weight;  
            % Covariance
            Ycov = K11-K12*model.net.invPrior*K12';
            yPredVar = Ycov;
            Kss_joint = K11;
        else
            [yPred, yPredVar] = gpfwd(model.net, X, model.cninv);
        end
        
    case 'BLR'
        % for Bayesian linear regression
        X = X(:, model.subs);
        [yPred, yPredVar] = bayes_fwd(model.mu, model.Sigma, model.beta, X);
    case 'LR'
        % for linear regression
        [yPred, yPredVar] = LRT(model, X);
        if model.options.logModel
            yPred = min(yPred, 99); % 10^99 is close to infinity in Matlab.
        end
    case 'smoother'
        yPred = zeros(size(X,1),1);
        for i=1:(size(X,1))
            distances = sq_dist(model.X',X(i,:)');
            weights = 1./(distances+0.001);
            weights = weights/sum(weights);
            yPred(i) = weights' * model.y;
        end
    case 'numOptW'
        X_prime = [ones(size(X,1),1), X];
        yPred = X_prime*model.w;
        yPredVar = zeros(size(yPred,1),size(yPred,2));
    case 'perturbInputs'
        yPredModels = zeros(size(X,1), model.numModels);
        for i=1:model.numModels
            yPredModels(:,i) = applyModel(model.submodel{i}, X);
        end
        yPred = mean(yPredModels,2);
        yPredVar = std(yPredModels,0,2).^2;
    case 'regression-tree'
        % How much to prune for prediction is a parameter to be tuned in cross-validation, params(2)
        params = (model.params - model.options.kernelParamLowerBound) / (model.options.kernelParamUpperBound-model.options.kernelParamLowerBound);
        numPruning = params(2) * (length(model.T.ntermnodes)-1);
        [yPred,nodes] = treeval(model.T, X, [length(model.T.ntermnodes) - numPruning]);
%        [yPred,nodes] = treeval(model.T, X, [length(model.T.ntermnodes)-10]);
%        [yPred,nodes] = treeval(model.T, X);
        yPredVar = ones(length(yPred),1);
        %yPredVar = model.T.nodeerr(nodes);
%             yPredVar = model.T.sigma(nodes).^2;
%             yExpImp = exp_imp_t(model.T.min_y, model.T.mu(nodes), model.T.sigma(nodes), model.T.nu(nodes));

    case {'randomregtree', 'regtree'}
        %yPred = fh_marg_treeval(model.T, X, model.options.transformation);
%        yPredPrime = treeval(model.T, X);
% more efficient if size(X,1) is large:       yPredPrime = fh_treeval_using_theta_pis(model.T, X);
         yPred = fh_simple_one_treeval(model.T, X, options.strategyForMissing);
%         assert(all(yPred<yPredPrime+0.0001) && all(yPred>yPredPrime-0.0001))
        yPredVar = ones(length(yPred),1);

        
    case 'ogp'
        [yPred, yPredVar] = ogpfwd(X);
%        stdT          = sqrt(varT);
%        global net;
%        meanBV        = ogpfwd(net.BV);

    otherwise
        if strfind(model.type, 'GP')
            global gprParams;
            gprParams = [];
            gprParams.combinedcat = model.cat;
            gprParams.combinedcont = model.cont;
%             gprParams.algoParam = model.algoParam;

            %=== Normalize X using same normalization as before:
%             X = X(:, model.good_feats);
%             X = X - repmat(model.means, [size(X,1),1]);
%             X = X ./ repmat(model.stds, [size(X,1),1]);


            if isfield(model.options, 'ppSize') && model.options.ppSize > 0
                if joint
                    error 'Joint predictions are not implemented yet for SRPP.'
                end
                    
                %=== GP with subset of regressors or projected process
                [yPred, S2SR, S2PP] = gprSRPPfwd(model.Kmm, model.invKmm, model.saved_1, model.saved_2, model.params, model.covfunc, model.pp_index, model.X, X, observationNoise);
                yPredVar = S2PP;
                if (any(yPredVar < 0))
                    debug_filename = 'debug_file_for_neg_var_in_pp.mat';
                    bout(sprintf(strcat(['\n\nWARNING: predicted variance is negative: ', num2str(min(yPredVar)), ', saving workspace to ', debug_filename])));
                    save(debug_filename);
                end
                yPredVar = max(yPredVar, 1e-10);
%                yPredVar = S2SR; % the two seem very similar, but the GP book suggests PP is better far away from the data.

            else
                %=== Normal GP
                if isfield(model, 'useCensoring')
                    if nargout == 3
                        [yPred, yPredVar] = gprCensorFwd(model.X, model.L_nonoise, model.alpha_nonoise, model.invK_times_invH_times_invK, model.covfunc, model.params, X, observationNoise, joint);
                    elseif nargout == 2
                        [yPred, yPredVar] = gprCensorFwd(model.X, model.L_nonoise, model.alpha_nonoise, model.invK_times_invH_times_invK, model.covfunc, model.params, X, observationNoise, joint);
                    else
                        yPred = gprCensorFwd(model.X, model.L_nonoise, model.alpha_nonoise, model.invK_times_invH_times_invK, model.covfunc, model.params, X, observationNoise, joint);
                    end
                 else
                    if nargout == 3
                        [yPred, yPredVar] = gprFwd(model.X, model.L, model.invL, model.alpha, model.covfunc, model.params, X, observationNoise, joint);

                        N_for_sampling_min = 200;
                        N_c = 100;
                        %=== Get joint prediction for N_for_sampling_min most promising
                        %design points; then take joint samples of those.
                        [tmp,idx] = sort(yPred-sqrt(yPredVar));
                        idx = idx(1:min(N_for_sampling_min,length(idx)));
                        [mean_best, tmp2, jointvar_best] = gprFwd(model.X, model.L, model.invL, model.alpha, model.covfunc, model.params, X(idx,:), observationNoise, 1);
                        samples = randnorm(N_c, mean_best, [], jointvar_best + eye(length(jointvar_best))*1e-8);
                        samples_of_min = min(samples,[],1)';
                        
                        %TODO: sample from joint of 50 best to get samples_of_min.
                    elseif nargout == 2
                        [yPred, yPredVar] = gprFwd(model.X, model.L, model.invL, model.alpha, model.covfunc, model.params, X, observationNoise, joint);
                    else
                        yPred = gprFwd(model.X, model.L, model.invL, model.alpha, model.covfunc, model.params, X, observationNoise, joint);
                    end
                end
            end

        else
            error ('No such model type defined yet!');
        end
end
if isfield(model, 'special_for_trivial') && (model.special_for_trivial==1)
    complete_yPred = log10(0.005) * ones(origSizeOfTest,1);
    complete_yPredVar = zeros(origSizeOfTest,1);
    complete_yPred(nontrivial_insts) = yPred;
    complete_yPredVar(nontrivial_insts) = yPredVar;
else
    complete_yPred = yPred;
    complete_yPredVar = yPredVar;
end