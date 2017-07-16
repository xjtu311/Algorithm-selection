function simple_prediction_experiment
% Carries out the experiment for simple empirical hardness models
% (i.e., only features vary).

k = 10; % k-fold crossval.


%options1 = get_default_activeconf_opts;
options_rf.modelType = 'rf';
options_rf.nSub = 10;
options_rf.logModel = 1;
options_rf.linearSize = 5;
options_rf.maxModelSize = 7;

options_rf.pca = 0; % do no PCA here, for feature importance
options_rf.init_split_ratio = 5.0/6; % options.tuning_params(3);
options_rf.init_Splitmin = 10; % options.tuning_params(2);
options_rf.orig_rf = 0;
options_rf.strategyForMissing = 0;
options_rf.doQuadratic = 1;

options_gp.modelType = 'GPML';
options_gp.pca = 0;
options_gp.logModel = 1;
options_gp.ppSize = 1000;
options_gp.trainSubSize = 1000;
options_gp.doQuadratic = 1;

options_lr.modelType = 'LR';
options_lr.pca = 0;
options_lr.logModel = 1;
options_lr.doQuadratic = 1;
options_lr.linearSize = 5;
options_lr.maxModelSize = 7;

options_vec = {options_lr,options_rf}%options_gp,options_rf,};
if isunix
    dir = '/ubc/cs/project/arrow/hutter/model_paper/';
else
    dir = 'Z:\arrowspace\model_paper\';
end

mip_benchmark_files = {};
sat_benchmark_files = {};
%mip_benchmark_files{end+1,1} = 'random-CPLEX-Orlib-results.txt';
% mip_benchmark_files{end+1,1} = 'scip-corlat.csv';
% mip_benchmark_files{end+1,1} = 'gurobi-corlat.csv';
% mip_benchmark_files{end+1,1} = 'cplex-corlat.csv';
% % mip_benchmark_files{end+1,1} = 'scip-bigmix.csv';
% % mip_benchmark_files{end+1,1} = 'gurobi-bigmix.csv';
mip_benchmark_files{end+1,1} = 'cplex-bigmix.csv';


% sat_benchmark_files{end+1,1} = '';

benchmark_files = [mip_benchmark_files; sat_benchmark_files];


for i=1:length(benchmark_files)
    benchmark_file = strcat(dir, benchmark_files{i});
    
    %=== Get X and y from that file.
    savedXy_file = strcat(benchmark_file, '-saved');
    if exist(savedXy_file, 'file')
        data = csvread(savedXy_file);
        X = data(:,1:end-1);
        y = data(:,end);
        featureNames = get_mip_feature_names;
    else
        y = csvread(benchmark_file,0,1); % runtimes

        instance_names = textread(benchmark_file,'%s%*[^\n]', 'bufsize', 10000);
        for j=1:length(instance_names)
            instance_names{j} = instance_names{j}(1:end-1); % remove the trailing comma
        end

        %=== Get features for instance names from database.
%         if mysql('status')
%             startUpDB
%         end
%         mysql('use MODEL_DB');
%         global feature_table;
%         feature_table = 'MIP_FEATURES';
%         inst_ids = get_inst_id(instance_names);
%         if i<= length(mip_benchmark_files)
%             mip = 1;
%         else
%             mip = 0;
%         end   
%         [X, featureNames] = getFeaturesForID(mip, inst_ids);
%         csvwrite(savedXy_file, [X, y]);

%         feature_file = strcat(dir, 'MIP_feats.csv');
% %        feature_file = strcat(dir, 'bigmix_feats.csv');
%         featureNames = {};
%         X = readFeatures(instance_names, feature_file);

    end
    
    %=== Shuffle data to avoid any biases in the cross-validation.
    perm = randperm(length(y));
%    X = X(:, [1:26,36:39]);

% TODO: get rid of this manual dropping of bad features!
    X = X(:,find(std(X)>0.0001));
    X = X(perm,:);
    y = y(perm);

    for j=1:length(options_vec)
        options = options_vec{j};
        
        %=== Do cross-validation to evaluate performance of different models.
        N = length(y);
        rand('twister',1234);

        startIdx = 1;
        totalLearnTime = 0;
        for l=1:k
            fprintf(strcat(['Cross-validation ', num2str(l), '/', num2str(k), '...\n']));
            endIdx = ceil(l*N/k);
            testIdx = startIdx:endIdx;
            trainIdx = setdiff(1:N, startIdx:endIdx);

            %trainIdx = trainIdx(1:300);
            trainIdx = trainIdx(1:1000);
            
            %=== Learn model for training data of this fold.
            cens = zeros(length(y),1);       % counting censored data at the censoring threshold.
            cat = [];                        % no categorical inputs
            catDomains = [];                 % still none
            tic;
            if options.logModel
                y = max(y, 0.005);
            end
            model = learnModel(X(trainIdx,:), y(trainIdx), cens(trainIdx), cat, catDomains, 0, options, featureNames);
            totalLearnTime = totalLearnTime + toc;
            
            %=== Predict for test data of this fold.
            [y_cross(testIdx,1), y_cross_var(testIdx,1)] = applyModel(model, X(testIdx,:), 0, 0, 0);

            startIdx = endIdx + 1;
        end

        if model.options.logModel
            [rmse, ll, cc] = measures_of_fit(log10(y), y_cross, y_cross_var, cens)
        else
            [rmse, ll, cc] = measures_of_fit(y, y_cross, y_cross_var, cens)
        end
        
        figure_prefix = strcat(dir, benchmark_files{i}, '-', options.modelType);
        title_prefix = '';
        plot_simple_pred_scatter(y, y_cross, y_cross_var, cens, rmse, cc, ll, figure_prefix, title_prefix, model.options.logModel)
        
        results{i,j} = [rmse, ll, cc, totalLearnTime];
    end    
end
results
for i=1:length(benchmark_files)
    fprintf(strcat(benchmark_files{i}, ' & '));
    two_output(results{i,1}(1),results{i,2}(1),1);
    fprintf(' & ');
%    two_output(results{i,1}(2),results{i,2}(2),0);
%    fprintf(strcat(benchmark_files{i}, ' & '));
    two_output(results{i,1}(3),results{i,2}(3),0);
    fprintf(' & ');
    two_output(results{i,1}(4),results{i,2}(4),1);
    fprintf(' \\\\\n');
end

function mip_feature_names = get_mip_feature_names
mip_feature_names = {'n_constraints'; 'n_variables'; 'nzcnt'; 'probtype'; 'num_q_constr'; 'num_quad'; 'num_qpnz'; 'vcg_var_deg_mean'; 'vcg_var_deg_std'; 'vcg_var_deg_min'; 'vcg_var_deg_max'; 'vcg_con_deg_mean'; 'vcg_con_deg_std'; 'vcg_con_deg_min'; 'vcg_con_deg_max'; 'prices_stddev'; 'price_per_number_stddev'; 'price_per_sqrt_number_stddev'; 'perc_cont_var'; 'support_size_mean'; 'support_size_std'; 'a_ij_magnitude_mean'; 'a_ij_magnitude_std'; 'a_ij_varcoef_avg'; 'a_ij_varcoef_std'; 'perc_unbounded_discrete'; 'edge_density'; 'vg_deg_std'; 'vg_deg_max'; 'vg_deg_min'; 'vg_deg_q0_25'; 'vg_deg_med'; 'vg_deg_q0_75'; 'clust_coef'; 'deviation'; 'lp_avg'; 'lp_l2_avg'; 'lp_linf'; 'lp_objval'};
%clean: [1:26,36:39]

function X = readFeatures(instance_names, feature_file)
X_file = csvread(feature_file, 0, 1);
instance_names_file = textread(feature_file,'%s%*[^\n]', 'bufsize', 200000);

for i=1:length(instance_names)
    for j=1:length(instance_names_file)
        succ=0;
        if strcmp(instance_names{i}, instance_names_file{j})
            X(i,:) = X_file(j,:);
            succ=1;
            break
        end
    end
    if(succ==0)
        instance_names{i}
        error(strcat(['no instance features not defined for instance ', instance_names{i}, ' in file ', feature_file]));
    end
end