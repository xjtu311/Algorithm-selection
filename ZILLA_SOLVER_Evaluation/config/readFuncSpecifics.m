function [func, rootdir] = readFuncSpecifics(scenario, numRun)
func.tuningScenario = scenario;

global feature_table;
% if strcmp(scenario, 'SPEAR-ibm-al') || strcmp(scenario, 'SPEAR-swv-al')
%     feature_table = 'FH_FEATURES_1sLS';
% elseif strcmp(scenario, 'SPEAR-competitions-al')% || strcmp(scenario, 'SAPS-SatComp-crafted-al') || strcmp(scenario, 'SAPS-SatComp-random-al')
%     feature_table = 'tmpFEATURES';
% %elseif strcmp(scenario, 'CPLEX-regions100') || strcmp(scenario, 'CPLEX-regions100-100train') || strcmp(scenario, 'CPLEX-regions200') 
% %    feature_table = 'NEW_FEATURES';
% elseif strcmp(scenario, 'CPLEX-cls-al') || strcmp(scenario, 'CPLEX-CLS-al') || strcmp(scenario, 'CPLEX-mik-al') || strcmp(scenario, 'CPLEX-MIK-al') || strcmp(scenario, 'CPLEX-conic-al') || strcmp(scenario, 'CPLEX-CONIC') || strcmp(scenario, 'CPLEX-conic5') || strcmp(scenario, 'CPLEX-regions100-al') || strcmp(scenario, 'CPLEX-regions100-100train-al') || strcmp(scenario, 'CPLEX-regions200-al') 
%     feature_table = 'MIP_FEATURES';
% else
%     feature_table = 'FEATURES'; %also for CPLEX instances from KLB's website
% end

%if strcmp(scenario, 'CPLEX-cls') || strcmp(scenario, 'CPLEX-CLS') || strcmp(scenario, 'CPLEX-mik') || strcmp(scenario, 'CPLEX-MIK') || strcmp(scenario, 'CPLEX-conic') || strcmp(scenario, 'CPLEX-CONIC') || strcmp(scenario, 'CPLEX-conic5') || strcmp(scenario, 'CPLEX-regions100') || strcmp(scenario, 'CPLEX-regions100_short') || strcmp(scenario, 'CPLEX-regions100-100train') || strcmp(scenario, 'CPLEX-regions200')  || strcmp(scenario, 'CPLEX-regions200_short') || strcmp(scenario, 'CPLEX-OR-al') || strcmp(scenario, 'CPLEX-cls-al_short') || strcmp(scenario, 'CPLEX-mik-a_short') || strcmp(scenario, 'CPLEX-conic-al_short') || strcmp(scenario, 'CPLEX-regions200-al_short') 
%     feature_table = 'MIP_FEATURES';
%else
%    feature_table = 'FEATURES'; %also for CPLEX instances from KLB's website
%end
% SAT table is always FH_THESIS_FEATURES, configured in getFeaturesForID()

if isempty(strfind(lower(scenario), 'cplex')) && isempty(strfind(lower(scenario), 'gurobi')) && isempty(strfind(lower(scenario), 'lpsolve')) && isempty(strfind(lower(scenario), 'lp_solve'))
    feature_table = 'NEWFEATURES';
    func.mip = 0;
else
    feature_table = 'MIP_FEATURES';
    func.mip = 1;
end



% global func
if nargin < 2
    numRun = 0;
end
func.runlength = 0;
func.scenario = scenario;
func.cutoff = 5;
func.tuningTime = 18000;
func.db = 1;
%func.overallobj = 'median';
%func.overallobj = 'mean10';
%func.overallobj = 'mean'; 
func.numTrainingInstances = 1000;
func.numTestInstances = 1000;

%=== Initialize with empty values, so don't assign if N/A.
func.cheap = 0;
func.matlab_fun = 0;
func.cat = [];
func.cont = [];
func.param_names = [];
func.default_values = [];
func.all_values = {};
func.param_bounds = [];

func.singleRunObjective = 'runtime';
func.overallobj = 'mean10'; % defaults unless otherwise specified in tuning scenario file.

unix_rootdir = '/.autofs/csother/ubccsprojectarrow/hutter/altuning/';
if isunix
    rootdir = unix_rootdir;
else
    rootdir = 'Z:\arrowspace\altuning\';
end
env.script_path = unix_rootdir;

%=== Read configuration from scenario file
filename = strcat(rootdir, 'tuning_scenarios/', scenario, '.txt');
if ~exist(filename, 'file')
    error(strcat(['Scenario file ', filename, ' does not exist.']));
end
lines=textread(filename, '%s', 'delimiter', '\n');
for i=1:length(lines)
    line = lines{i};
    matches = regexp(line, '(.*)=(.*)', 'tokens', 'once');
    if length(matches) ~= 2
        error(strcat('Line ', line, ' does not match syntax <param>=<value>'));
    end
    name = ddewhite(matches{1});
    value = ddewhite(matches{2});
    
    switch name
        case 'algo'
            env.algo = value;
%             if strcmp(env.algo, 'saps')
%                 env.algo_wrapper = 'al_saps_wrapper.rb'
%             elseif strcmp(env.algo, 'spear')
%                 env.algo_wrapper = 'al_spear_wrapper.rb';
%             elseif strcmp(env.algo, 'cplex')
%                 error 'Must figure out the wrapper for CPLEX first'
%             else
%                 error 'Unknown algorithm -- change readFuncSpecifics.rb and all scenariofiles to include env.algo_wrapper for this algorithm'
%             end
        case 'exec_path'
            env.exec_path = strcat([unix_rootdir, value]);
        case 'params_filename'
            func.params_filename = strcat(rootdir, value);
            env.params_filename = strcat([unix_rootdir, value]);
        case 'deterministic'
            func.deterministic = str2num(value);
        case 'cheap'
            func.cheap = str2num(value);
        case 'matlab_fun'
            func.matlab_fun = str2num(value);
        case 'executable'
            env.executable = value; %strcat(unix_rootdir, value);
        case 'instance_seed_file_prefix'
            local_instance_seed_file_prefix = strcat(rootdir, value);
            remote_instance_seed_file_prefix = strcat(unix_rootdir, value);
        case 'test_instance_seed_file'
            local_test_instance_seed_file = strcat(rootdir, value);
            remote_test_instance_seed_file = strcat(unix_rootdir, value);
        case 'outdir'
            if isunix
                delim = '/';
            else
                delim = '\';
            end
            func.outdir = strcat([rootdir, value, delim]);
            func.remote_outdir = strcat([unix_rootdir, value, '/']);
        case 'tunerTimeout'
            func.tuningTime = str2num(value);
        case 'cutoff'
            func.cutoff = str2num(value);
        case 'numTrainingInstances'
            func.numTrainingInstances = str2num(value);
        case 'numTestInstances'
            func.numTestInstances = str2num(value);
        case 'runlength'
            func.runlength = str2num(value);
        case 'db'
            func.db = str2num(value);
        case 'singleRunObj'
            func.singleRunObjective = value;
        case 'overallobj'
            func.overallobj = value;
            
            
        otherwise
            error(strcat('Unknown configuration parameter ', name, ' in line ', line))
    end
end
func.env = env;
func.outdirScen = strcat([func.outdir, scenario, '/']);

% %=== For deployment, can split this in one file for each scenario
% switch scenario
%     case 'SWa-SAPS'
%         func.name = 'saps'
%         func.transformation = 'log';
%         func.cheap = 0;
%
%         %=== Same for all SAPS scenarios.
%         env.algo = 'saps';
%         env.executable = 'ruby saps_wrapper.rb';
%         func.params_filename = 'tuning_scenarios/saps/saps-params.txt';
%
%     case 'SWa-SAPS-cont'
%         func.name = 'saps'
%         func.transformation = 'log';
%         func.cheap = 0;
%
%         %=== Same for all SAPS scenarios.
%         env.algo = 'saps';
%         env.executable = 'ruby saps_wrapper.rb';
%
%         func.params_filename = 'tuning_scenarios/saps/saps-cont-params.txt';
%
%     case 'QCP-SAPS-cont'
%         func.name = 'saps'
%         func.transformation = 'log';
%         func.cheap = 0;
%
%         %=== Same for all SAPS scenarios.
%         env.algo = 'saps';
%         env.executable = 'ruby saps_wrapper.rb';
%
%         func.params_filename = 'tuning_scenarios/saps/saps-cont-params.txt';
%
%     case 'SWa-SPEAR-cont'
%         func.name = 'spear'
%         func.transformation = 'log';
%         func.cheap = 0;
%
%         %=== Same for all SAPS scenarios.
%         env.algo = 'spear';
%         env.executable = 'ruby spear_wrapper.rb';
%
%         func.params_filename = 'tuning_scenarios/spear/spear-cont-params.txt';
%
%     case 'QCP-SPEAR-cont'
%         func.name = 'spear'
%         func.transformation = 'log';
%         func.cheap = 0;
%
%         %=== Same for all SAPS scenarios.
%         env.algo = 'spear';
%         env.executable = 'ruby spear_wrapper.rb';
%
%         func.params_filename = 'tuning_scenarios/spear/spear-cont-params.txt';
%
%     case 'CPLEX-CATS-cont'
%         func.name = 'cplex'
%         func.transformation = 'log';
%         func.cheap = 0;
%
%         %=== Same for all CPLEX scenarios.
%         env.algo = 'cplex';
%         env.executable = 'ruby cplex_wrapper.rb';
%
%         func.params_filename = 'tuning_scenarios/cplex/cplex-cont-params.txt';
%
%     case 'SWa-SAPS-0-surrogate'
% %        load surrogateModelForSW_SAPS_0;
% %        load mixtureOfRegtreeSurrogateModelForSW_SAPS_0;
%         load smoothSurrogateModelForSW_SAPS_0;
%         func.name = 'saps';
%         func.transformation = 'log';
%         func.cheap = 1;
%         func.f_lb = 0;
%
%         %       func.funcHandle = @(x) applyModel(surrogateModelForSW_SAPS_0, x)
%         func.funcHandle = @(x) applyModel(model, x);
%         func.params_filename = 'tuning_scenarios/saps/saps-cont-params.txt';
%
%     case 'hartman6_marginal'
%         func.name = 'hartman';
%         func.cheap = 1;
%         func.funcHandle = @hartman6_marginalized;
%         func.transformation = 'id';
%
%         func.cont = 1:4;
%         func.cont_param_names = {'1','2','4','6'};
%         func.cont_param_bounds = [0,1;0,1;0,1;0,1];
%         func.cont_default_values = 1/2*(func.cont_param_bounds(:,2)+func.cont_param_bounds(:,1));
%         func.best_values = [0.38928; 0.87683; 0.58822; 0.03835]; %-1.13630;
%
% %     case 'hart6'
% %         func.name = 'hart6'
% %         func.cheap = 1;
% %         func.funcHandle = @hart6;
% %         func.dim = 6;
% %         func.transformation = 'id';
% %         func.lower = zeros(func.dim,1);
% %         func.upper = ones(func.dim,1);
% %         func.cat = [];
% %         func.sizes = zeros(func.dim,1);
% %         func.values = zeros(func.dim,1);
% %         func.best_values = [0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573];
%
% %     case 'hartman6_marginal_t'
% %         func.name = 'hartman'
% %         func.cheap = 1;
% %         func.funcHandle = @hartman6_marginalized_t;
% %         func.dim = 4;
% %         func.transformation = 'id';
% %         func.lower = ones(func.dim,1);
% %         func.upper = 2*ones(func.dim,1);
% %         func.cat = [];
% %         func.sizes = zeros(func.dim,1);
% %         func.values = zeros(func.dim,1);
% %         func.best_values = [0,0,0.58822,0.03835];
%
%     case 'braninprod'
%         func.name = 'braninprod'
%         func.cheap = 1;
%         func.funcHandle = @true_prod_of_branins;
%         func.transformation = 'log';
%
%         func.cont = 1:2;
%         func.cont_param_names = {'x1','x4'};
%         func.cont_param_bounds = [0,1;0,1];
%         func.cont_default_values = 1/2*(func.cont_param_bounds(:,2)+func.cont_param_bounds(:,1));
%         func.best_values = [0.20263; 0.25445];

% %      case 'gold'
% %         func.name = 'gold'
% %         func.cheap = 1;
% %         func.funcHandle = @gold;
% %         func.dim = 2;
% %         func.transformation = 'id';
% %         func.lower = [-2,-2]';
% %         func.upper = [2,2]';
% %         func.cat = [];
% %         func.sizes = zeros(func.dim,1);
% %         func.values = zeros(func.dim,1);
%
%     case 'sixHumpCamelBack'
%         func.name = 'sixHumpCamelBack';
%         func.cheap = 1;
%         func.funcHandle = @sixHumpCamelBack;
%         func.transformation = 'log';
%
%         func.cont = 1:2;
%         func.cont_param_names = {'x','y'};
%         func.cont_param_bounds = [-3,3;-2,2];
%         func.cont_default_values = 1/2*(func.cont_param_bounds(:,2)+func.cont_param_bounds(:,1));
%         func.best_values = [-0.0898; 0.7126];
%
% %     case 'branin_alg'
% %         func.name = 'branin_alg'
% %         func.cheap = 0;
% %         func.funcHandle = @branin_alg_fixed;
% %         func.dim = 2;
% %         func.transformation = 'id';
% %         func.lower = zeros(func.dim,1);
% %         func.upper = ones(func.dim,1);
% %         func.cat = [];
% %         func.sizes = zeros(func.dim,1);
% %         func.values = zeros(func.dim,1);
%
%     case 'disc_test_fun'
%         func.name = 'disc_test_fun'
%         func.cheap = 1;
%         func.funcHandle = @disc_test_fun;
%         func.dim = 4;
%         func.transformation = 'log';
%         func.cat = 1:func.dim;
%         func.sizes = 7*ones(func.dim,1);
%         for i=1:func.dim
% %             func.cat_param_names{i} = num2str(i);
%             func.cat_all_values{i} = 1:7;
%             func.cat_default_values{i} = 4;
%         end
%
% %     case 'hybrid_test_fun'
% %         func.name = 'simplest_discrete'
% %         func.cheap = 1;
% %         func.funcHandle = @simplest_discrete;
% %         func.dim = 4;
% %         func.transformation = 'id';
% %         func.cat = 1:2;
% %         func.sizes = 7*ones(func.dim,1);
% %         func.values = {};
% %         func.values{1} = 1:7;
% %         func.values{2} = 1:7;
% %         func.values{3} = 1:7;
% %         func.values{4} = 1:7;
% %         for i=1:func.dim
% %             func.lower(i) = min(func.values{i});
% %             func.upper(i) = max(func.values{i});
% %         end
%
%     case 'cont_test_fun'
%         func.name = 'simplest_discrete'
%         func.cheap = 1;
%         func.funcHandle = @simplest_discrete;
%         func.transformation = 'log';
%
%         func.cont = 1:10;
%         for i=1:length(func.cont)
%             func.cont_param_names{i} = num2str(i);
%             func.cont_param_bounds(i,:) = [0, 2];
%             func.cont_default_values = 1/2*(func.cont_param_bounds(:,2)+func.cont_param_bounds(:,1));
%             func.best_values = 2*ones(length(func.cont),1);
%             func.best_values(i) = 1;
%         end
%
%     case 'disc_test_fun2'
%         func.name = 'disc_test_fun2'
%         func.cheap = 1;
%         func.funcHandle = @disc_test_fun2;
%         func.transformation = 'log';
%
%         dim = 4;
%         func.cat = 1:dim;
%         for i=1:dim
% %             func.cat_param_names{i} = num2str(i);
%             func.cat_all_values{i} = 1:7;
%             func.cat_default_values{i} = 4;
%         end
%         % opt = 3.5795;
%
%     case 'disc_bnet_ll'
%         func.name = 'disc_bnet_ll'
%         func.cheap = 1;
%         func.funcHandle = @disc_bnet_ll;
%         bnet = mk_alarm_bnet;
%         func.transformation = 'id';
%
%         func.cat = 1:length(bnet.node_sizes);
%         for i=1:length(bnet.node_sizes)
% %             func.cat_param_names{i} = num2str(i);
%             vals = {};
%             for j=1:bnet.node_sizes(i)
%                  vals{j} = j;
%             end;
%             func.cat_all_values{i} = vals;
%             func.cat_default_values{i} = 2;
%         end
%         best = [1 1 1 1 1 2 2 1 1 1 1 1 1 2 3 2 1 2 2 2 3 2 3 3 3 3 2 2 3 2 1 2 2 2 2 2 2]';
%         func.best_values = {};
%         for i =1:length(best)
%             func.best_values{i} = best(i);
%         end
%
%     case 'saved_saps_swgcp_first100train'
%         func.name = 'saved_saps_swgcp_first100train';
%         func.cheap = 1;
%         func.transformation = 'log';
%
%         func.cat = 1:4;
% %         func.cat_param_names = {'alpha','rho','ps','wp'};
%         func.cat_all_values{1} = [1.01, 1.066, 1.126, 1.189, 1.256, 1.326, 1.4];
%         func.cat_all_values{2} = [0, 0.17, 0.333, 0.5, 0.666, 0.83, 1];
%         func.cat_all_values{3} = [0, 0.033, 0.066, 0.1, 0.133, 0.166, 0.2];
%         func.cat_all_values{4} = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06];
%         func.cat_default_values = [4;4;4;4];
%
%         func.best_values = [5,2,1,5]; %        [1.256, 0.17, 0, 0.04]
%         %=== Set up anonymous function.
%         func.funcHandle = @(x) saved_saps_swgcp_first100train(x, func.cat_all_values);
%
%     case 'simplest_discrete'
%         func.name = 'simplest_discrete'
%         func.cheap = 1;
%         func.funcHandle = @simplest_discrete;
%         func.transformation = 'log';
%
%         dim = 30;
%         func.cat = 1:dim;
%         for i=1:dim
% %             func.cat_param_names{i} = num2str(i);
%             func.cat_all_values{i} = {'0', '1', '2'};
%             func.cat_default_values{i} = '1';
%             func.best_values{i} = 3;
%             func.best_values{1} = 2;
%         end
%
% %     case 'mixed_disc_cont_testfun'
% %         func.name = 'mixed_disc_cont_testfun'
% %         func.cheap = 1;
% %         func.funcHandle = @mixed_disc_cont_testfun;
% %         func.dim = 4;
% %         func.transformation = 'id';
% %         func.cat = 1:2;
% %         func.sizes = 2*ones(func.dim,1);
% %         func.values = {};
% %         for i=1:func.dim
% %             func.values{i} = [-0.5,1];
% %             func.lower(i) = min(func.values{i});
% %             func.upper(i) = max(func.values{i});
% %         end
%
%     otherwise
%         error(strcat(['Unknown scenario: ' scenario]));
% end

%=== Same for all scenarios.
if ~exist(func.outdir, 'dir')
    mkdir(func.outdir);
end
if ~exist(func.outdirScen, 'dir')
    mkdir(func.outdirScen);    
end

%func.outputfilenameprefix = strcat([func.outdir, 'B-algo', env.algo, '-runobjruntime-overallobj', func.overallobj, '-runs', num2str(func.N), '-time', num2str(func.cutoff)]);
func.outputfilenameprefix = strcat([func.outdir, 'P', num2str(func.cutoff)]);

func.seeds = cell(func.numTrainingInstances,1);
if func.cheap || func.matlab_fun
    %=== Optimizing some function.
    func.outputfilenameprefix = strcat(func.outdir, 'B-time', num2str(func.cutoff));
    func.seeds = ceil((2^32-1)*rand(func.numTrainingInstances,100000));
    func.test_seeds = ceil((2^32-1)*rand(func.numTestInstances,100000));
else
    %=== Optimizing an algorithm -- read in instances, seeds, and info on the parameters.
    local_instance_seed_filename = strcat([local_instance_seed_file_prefix num2str(numRun) '.txt']);
    remote_instance_seed_filename = strcat([remote_instance_seed_file_prefix num2str(numRun) '.txt']);

    if func.deterministic
        [func.seeds, func.instance_filenames] = read_seed_instance_file(local_instance_seed_filename);
    else
        [func.seeds, func.instance_filenames] = read_instances_and_seeds_rnd(local_instance_seed_filename);
    end
% for debug    func.seeds = func.seeds(:, 1);
    
%    seeds = seeds(1:func.N);
%    instance_filenames = instance_filenames(1:func.N);
    
    [func.test_seeds, func.test_instance_filenames] = read_seed_instance_file(local_test_instance_seed_file);
    if length(func.test_seeds) > 1000
        func.test_seeds = func.test_seeds(1:1000);
        func.test_instance_filenames = func.test_instance_filenames(1:1000);
    end
    func.numTestInstances = length(func.test_instance_filenames);
%    test_seeds = seeds(1:func.N);
%    instance_filenames = instance_filenames(1:func.N);
end

if isfield(func,'params_filename')
    [func.cat, func.cont, func.param_names, func.all_values, func.param_bounds, func.param_trafo, func.is_integer_param, func.default_values, func.cond_params_idxs, func.parent_param_idxs, func.ok_parent_value_idxs] = read_params(func.params_filename);
    func.dim = length(func.cat)+length(func.cont);
    func.param_names = func.param_names';
    func.orig_param_lower_bound = func.param_bounds(:,1);
    func.orig_param_upper_bound = func.param_bounds(:,2);
    func.param_bounds(find(func.param_bounds(:,1)),1) = 0;
    func.param_bounds(find(func.param_bounds(:,2)),2) = 1;
    func.default_values = config_transform(func.default_values', func)';
else
    error('Need to specify params_filename in scenario file.')
end

%=== Get number of values for each parameter.
func.num_values = zeros(1,func.dim);
for i=1:func.dim
    func.num_values(i) = length(func.all_values{i});
end

%=== Compute number of free parameters.
func.effdim = 0;
for i=1:length(func.cont)
    j = func.cont(i);
    if func.param_bounds(j,1) < func.param_bounds(j,2)
        func.effdim = func.effdim+1;
    end
end
for i=1:length(func.cat)
    j = func.cat(i);
    if func.num_values(j)>1
        func.effdim = func.effdim+1;
    end
end

%=== If func.best_values not defined, use the defaults.
if ~isfield(func, 'best_values')
    func.best_values = func.default_values;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING UP ANONYMOUS FUNCTIONS - HAS TO BE IN THE END since func is passed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if func.cheap || func.matlab_fun
    %=== Anonymous function to evaluate the function read in as a string.
    handle = @(Theta, instance_numbers, seeds, censorLimits) feval(env.algo, Theta, instance_numbers, seeds, censorLimits);

    %=== Build anonymous function dealing with censoring.
    func.funcHandle = @(Theta_idx,instanceNumbers,seeds,censorLimits,name) censoredReading(func, handle, Theta_idx, instanceNumbers, seeds, censorLimits, func.cutoff, name, 0); % add noise!  1
    func.testFuncHandle = func.funcHandle;
%    func.noisefreeHandle = @(Theta_idx,instance_numbers,censorLimits,name) censoredReading(handle, Theta_idx, instance_numbers, func.seeds, censorLimits, func.cutoff, name, 0);
else
    %=== Set up anonymous function with complete struct func (thus at the end).
    %    func.funcHandle = @(x) eval_algo(func, x, instance_filenames, seeds, func.cutoff, func.db);
    %    func.funcHandle = @(x) get_batch_results(func, x, instance_filenames, seeds, instance_seed_filename);
    %    func.funcHandle = @(x,censorLimit) eval_algo(func, x, instance_filenames, seeds, instance_seed_filename, censorLimit);

    %    func.funcHandle = @(Theta,instance_filenames,seeds,censorLimits) get_single_results(func, Theta, instance_filenames, seeds, censorLimits);
    func.funcHandle = @(Theta_idx,instance_numbers,seeds,censorLimits,name) get_single_results(func, Theta_idx, instance_numbers, seeds, func.instance_filenames(instance_numbers), censorLimits, func.cutoff, func.runlength, name);
    func.testFuncHandle = @(Theta_idx,instance_numbers,seeds,censorLimits,name) get_single_results(func, Theta_idx, instance_numbers, seeds, func.test_instance_filenames(instance_numbers), censorLimits, func.cutoff, func.runlength, name);
end



function [reading, censored] = censoredReading(func, funcHandle, Theta_idx, instance_numbers, seeds, censorLimits, kappa, name, observationNoise)
global ThetaUniqSoFar;
global TestTheta;
global resfid;
global runTimeForRunningAlgo;

reading = -ones(length(Theta_idx),1);
censored = -ones(length(Theta_idx),1);
for i=1:length(Theta_idx)
    theta_idx = Theta_idx(i);
    if theta_idx < 0
        theta = TestTheta(-theta_idx,:);
    else
        theta = ThetaUniqSoFar(theta_idx,:);
    end
    theta = config_back_transform(theta, func);
    uncensoredReading = funcHandle(theta, instance_numbers(i), seeds(i), censorLimits(i));
    uncensoredReading = uncensoredReading;% + observationNoise*randn;
    if uncensoredReading < censorLimits(i)
        reading(i) = uncensoredReading;
        censored(i) = 0;
    else
        reading(i) = censorLimits(i);
        censored(i) = 1;
    end
    runTimeForRunningAlgo = runTimeForRunningAlgo + reading(i);
    fprintf(resfid, '%g, %d, %d, %d, %d', [reading(i), censored(i), theta_idx, instance_numbers(i), seeds(i)]);
    
    for j=1:size(theta,2)
        fprintf(resfid, ', %g', theta);
    end
    fprintf(resfid, '\n');
end
%reading = uncensoredReading; censored=zeros(size(X,1),1);