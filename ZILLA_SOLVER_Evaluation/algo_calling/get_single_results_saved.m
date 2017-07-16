function [ys, censoreds, runtimes, runlengths, solveds, best_sols] = get_single_results(func, Theta_idx, instance_numbers, seeds, instance_filenames, censorLimits, kappa, dorunlength, name, numrun)
%=== Returns the result (runtime/runlength) of the algorithm for the parameter configurations indexed by Theta_idx on the given instances with the given seeds.
%=== If using the database, get the appropriate IDs and check if the run is
%=== already stored there.
%=== Otherwise, call a script via the command line that executes the runs.

global TestTheta;
global ThetaUniqSoFar;

if nargin < 9
    error 'need 9 arguments' 
    %name = 'manyruns';
end
ys = -ones(length(Theta_idx),1);
censoreds = -ones(length(Theta_idx),1);
runtimes = -ones(length(Theta_idx),1);
solveds = -ones(length(Theta_idx),1);
runlengths = -ones(length(Theta_idx),1);
best_sols = -ones(length(Theta_idx),1);

oldLeftTodo = 1:length(Theta_idx);
assert(~isempty(oldLeftTodo));
numTry = 0;
Theta_idx_to_run = [];
insts_to_run = {};
seeds_to_run = [];
censorLimits_to_run = [];
files_to_delete = {};
cpu_times = -1e10*ones(length(Theta_idx),1);
while ~isempty(oldLeftTodo) && numTry<20
    numTry = numTry+1;
    if numTry > 2
        numTry=numTry
    end
    leftTodo = oldLeftTodo;
    for j = 1:length(oldLeftTodo)
        if (mod(j, 1000)==0)
            j
        end
        i = oldLeftTodo(j);
        theta_idx = Theta_idx(i);
        instance_filename = instance_filenames{i};
        seed = seeds(i);
        censorLimit = censorLimits(i);
        if func.db
            if any(theta_idx<0) % validation configs.
                Theta = TestTheta(-theta_idx,:);
            else
                Theta = ThetaUniqSoFar(theta_idx,:);
            end
            [runtime, censored, runlength, best_sol, solution_status] = ask_db(func, Theta, instance_filename, censorLimit, seed, dorunlength);
        else
            [runtime, censored, runlength, best_sol, solution_status, files_to_delete{end+1}] = doRunAndGetResults(func, theta_idx, instance_filename, censorLimit, seed, dorunlength);
        end
        
        if runtime > 10
            runtime
        end
        
        
        if dorunlength
            cpu_time = min(censorLimit, runlength);
        else
            cpu_time = min(censorLimit, runtime);
        end

        if isempty(runtime)
            if length(oldLeftTodo) == 1
                if dorunlength
                    startRun(func, theta_idx, instance_filename, seed, 5, 0);
                else
                    startRun(func, theta_idx, instance_filename, seed, censorLimit, 0);
                end
            else
                i
                %=== Only wait for evaluation of last point.
    %             waitToFinish = 0;
    %             if i == max(oldLeftTodo)
    %                 waitToFinish = 1;
    %             end
                % startRun(func, theta_idx, instance_filename, seed, censorLimit, waitToFinish) 
                Theta_idx_to_run(end+1,:) = theta_idx;
                insts_to_run{end+1} = instance_filename;
                seeds_to_run(end+1) = seed;
                if dorunlength
                    censorLimits_to_run(end+1) = 5;
                else
                    censorLimits_to_run(end+1) = censorLimit;
                end
            end
        else
            runlengths(i) = runlength;
            runtimes(i) = runtime;
            censoreds(i) = censored;
            cpu_times(i) = cpu_time;
            solveds(i) = solution_status;
            best_sols(i) = best_sol;

            leftTodo = setdiff(leftTodo, i);
        end
        if ~isempty(insts_to_run) && mod(length(insts_to_run), 10001)==0 && j < length(oldLeftTodo)
            % Start runs, don't wait and empty queque of runs to do.
            files_to_delete{end+1} = startManyRuns(func, Theta_idx_to_run, insts_to_run, seeds_to_run, censorLimits_to_run, 1, name);
            Theta_idx_to_run = [];
            insts_to_run = {};
            seeds_to_run = [];
            censorLimits_to_run = [];
        end
    end
    if ~isempty(insts_to_run)
        files_to_delete{end+1} = startManyRuns(func, Theta_idx_to_run, insts_to_run, seeds_to_run, censorLimits_to_run, 1, strcat('L-', name)); 
        Theta_idx_to_run = [];
        insts_to_run = {};
        seeds_to_run = [];
        censorLimits_to_run = [];
    end
	oldLeftTodo = leftTodo
end
if ~isempty(oldLeftTodo)
    func
    Theta_idx
    instance_filenames
    seeds
    censorLimits
    name
    oldLeftTodo
    numTry
    insts_to_run
    error('Still not successful after 20 tries.');
end
global runTimeForRunningAlgo;
assert(all(cpu_times >= -1e-5));
runTimeForRunningAlgo = runTimeForRunningAlgo + sum(cpu_times);

%=== Clean up files in directory RUNFILES - can't do earlier due to weird
%=== error (probably since Matlab is looking ahead, executing the delete
%=== too early?)
for i=1:length(files_to_delete)
    fprintf(strcat('Deleting file: ', files_to_delete{i}));
    delete(files_to_delete{i});
end

switch func.singleRunObjective
    case 'runtime'
        ys = runtimes;
    case 'runlength'    
        ys = runlengths;
    case 'solqual'
        ys = best_sols;
        censoreds = zeros(size(y,1),1);
    otherwise 
        error 'Still have to implement objective functions for single algorithm runs other than runtime, runlength, and solqual.'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% HELPER FUNCTIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORIGINAL VERSION BEFORE CHANGING TO WRITING INTO DB FROM WITHIN MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function startRun(func, theta_idx, instance_filename, seed, censorLimit, waitToFinish)
env = func.env;
%=== Construct param_string.
global allParamStrings;
param_string = allParamStrings{theta_idx};
%alphabeticalParameterString(func, theta);

%=== Call run via the command line.
run_cmd = strcat(['ruby scripts/single_runstarter.rb ' env.algo ' ' env.exec_path ' ' env.params_filename ' ' instance_filename ' ' num2str(seed) ' ' num2str(censorLimit) ' 0 1 ' param_string]);
if isunix
    connect_cmd = strcat(['cd ' env.script_path '; ']);
    cmd = strcat([connect_cmd ' ' run_cmd])
    %[a,b] = unix(cmd);
    call(cmd); % faster version, but no outputs
else
    %                error('Dont want to run algos from DOS')
    connect_cmd = strcat(['ruby remote_executer.rb ' env.script_path]);
    cmd = strcat([connect_cmd ' ' quote(run_cmd)]);
    [a,b] = dos(cmd);
end


function [runtime, censored, runlength, best_sol, solution_status, remote_file_to_delete] = doRunAndGetResults(func, theta_idx, instance_filename, censorLimit, seed, dorunlength)
global numrun
if theta_idx < 0
    global TestTheta;
    param_string = alphabeticalParameterString(func, TestTheta(-theta_idx,:));
else
    global allParamStrings;
    param_string = allParamStrings{theta_idx};
end
env = func.env;

%=== Replicating the hack from get_runtime.m
if dorunlength
    cutoff_length = censorLimit; 
    censorLimit = 5;% HACK since we only have an interface to specify cutoff_time not cutoff_length -> in this time we have to do more steps than cutoff_length.
else
    cutoff_length = 1000000;
end

% ruby saps_wrapper.rb <instance_relname> <instance_specifics> <cutoff_time> <cutoff_length> <seed> <params to be passed on>."

%=== Call run via the command line.
if isunix
    s = rand('state');
    deldir = strcat(func.outdir, 'DEL2RUNFILES/');
    subdir = strcat(deldir, num2str(numrun), '/');
    if ~exist(deldir, 'dir')
        mkdir(deldir);
    end
    if ~exist(subdir)
        mkdir(subdir);
    end

    outfiledir = subdir;
    %outfile = strcat(env.exec_path, '/', env.algo, random_string(), '-out.txt');
    outfile = strcat(outfiledir, env.algo, random_string(), '-out.txt');
    run_cmd = strcat([env.executable ' ' instance_filename ' 0 ' num2str(censorLimit) ' max ' num2str(seed) ' ' outfile ' ' param_string]);
    connect_cmd = strcat(['cd ' env.exec_path '; ']);
    cmd = strcat([connect_cmd ' ' run_cmd])
    status = 1;
    numtries = 1;
    while status ~= 0 && numtries <= 20
%         [status, result] = unix(cmd);
        call(cmd);
        status = ~exist(outfile, 'file'); % 0 if the outfile exists.
        if status ~= 0
            bout(strcat(['Try number ', num2str(numtries), ' failed for running command ', cmd]));
        end
        numtries = numtries + 1;
    end
    %=== Just parse result and return -> not writing to DB here.
%     res = regexp(result, ',', 'split');

    %=== Not parsing output directly, but letting wrapper output the relevant info to a file.
    res = csvread(outfile);
    solution_status = res(end-4);
%     switch solution_status
%         case 0
%             solution_status = 'TIMEOUT';
%         case 1
%             solution_status = 'SAT';
%         case 2
%             solution_status = 'UNSAT';
%         case 3 
%             solution_status = 'WRONG ANSWER';
%         otherwise
%             error 'unknown solution status'
%     end
    if solution_status == 3 % wrong answer!
        warnstring = strcat(['Wrong answer for theta_idx ', num2str(theta_idx), ', instance ', num2str(instance_filename), ', and seed ', num2str(seed), '.\nBreaking command:', cmd]);
        warning(warnstring);
        bout(warnstring);
        
        seed_out = res(end);    
        assert(seed == seed_out);
        runtime = inf;
        runlength = inf;
        best_sol = inf;
    else
        runtime = res(end-3);
        runlength = res(end-2);
        best_sol = res(end-1);
        seed_out = res(end);    
        assert(seed == seed_out);
    end
    
    fprintf(strcat('Deleting file: ', outfile));
    remote_file_to_delete = outfile;
    
    if strcmp(solution_status, 'TIMEOUT')
        assert(runtime >= censorLimit-1e-6);
        censored = 1;
        runtime = censorLimit;
        if dorunlength
            runlength = cutoff_length;
        end
        return
    else
        censored = 0;
    end
else
    error 'Running algorithms directly from DOS not implemented -> have to use the DB or run from UNIX.'
end

function remote_file_to_delete = startManyRuns(func, Theta_idx, instance_filenames, seeds, censorLimits, waitToFinish, name) 
if nargin < 7
    name = 'manyruns';
end

% Don't want to rewrite all the functionality of adding inst_id,
% algo_config_id (!) and algorun_config_id to the database.
% Rather use the existing ruby scripts to do that, we only use
% Matlab for queries: simply write a file where each line is a
% call to single_runstarter.rb
env = func.env;
global allParamStrings;
global clusterize;
global numrun
subRunClusterize = 0;
run_cmds = {};

%=== Create file with commands to run.
s = rand('state');
deldir_post = strcat('DEL3RUNFILES/');
subdir_post = strcat(deldir_post, num2str(numrun), '/');
if ~exist(strcat(func.outdir, '/', deldir_post), 'dir')
    mkdir(strcat(func.outdir, '/', deldir_post));
end
if ~exist(strcat(func.outdir, '/', subdir_post), 'dir')
    mkdir(strcat(func.outdir, '/', subdir_post));
end

file_with_cmds = strcat(subdir_post,strrep(strrep(datestr(now), ' ', '-'), ':', '-'), '--', num2str(rand));
remote_file_with_cmds = strcat(func.remote_outdir, file_with_cmds);
file_with_cmds = strcat(func.outdir, file_with_cmds);
%if ~exist(strcat(func.outdir, 'DELRUNFILES/'), 'dir')
%    mkdir(strcat(func.outdir, 'DELRUNFILES/'));
%end
rand('state', s);


%=== Run the commands in the file, if clusterize then on the cluster.
if clusterize
    if waitToFinish>0
        oncluster = 6;
    else
        oncluster = 8;
    end
else
    oncluster = 0;
end

global TestTheta;
cmds_fid = fopen(file_with_cmds,'w');
fprintf(cmds_fid, strcat([env.algo, ' ', env.exec_path, ' ', env.params_filename, ' ', num2str(oncluster), '\n']));
for i=1:length(Theta_idx)
    %=== Construct param_string.
    theta_idx = Theta_idx(i);
    if theta_idx < 0
        param_string = alphabeticalParameterString(func, TestTheta(-theta_idx,:));
    else
        param_string = allParamStrings{theta_idx};
    end
    fprintf(cmds_fid, strcat([instance_filenames{i}, ' ', num2str(num2str(seeds(i))), ' ', num2str(censorLimits(i)), ' ', param_string, '\n']));
end
fclose(cmds_fid); % throws error if we're deleting the file after.
run_cmd = strcat(['ruby scripts/al_run_configs_in_file.rb ', remote_file_with_cmds]) 
remote_file_to_delete = remote_file_with_cmds;

if isunix
    connect_cmd = strcat(['cd ' env.script_path '; ']);
    cmd = strcat([connect_cmd ' ' run_cmd]);
    [a,b] = unix(cmd);
    a
    b
else
    %                error('Dont want to run algos from DOS')
    connect_cmd = strcat(['ruby remote_executer.rb ' env.script_path]);
    cmd = strcat([connect_cmd ' ' quote(run_cmd)])
    [a,b] = dos(cmd);
end
%delete(file_with_cmds); % makes fclose above throw an error occasionally


%=== Get runtime from DB (or runlength, if input runlength == 1).
function [runtime, censored, runlength, best_sol, solution_status] = ask_db(func, theta, inst_name, cutoff, seed, dorunlength)
runtime = [];
censored = [];
runlength = [];
best_sol = [];
solution_status = [];
possible_inst_ids = get_inst_id({inst_name});
possible_inst_ids = possible_inst_ids{1};
for i=1:length(possible_inst_ids)
    inst_id = possible_inst_ids(i);
    if ~isempty(inst_id)
        algo_config_id = get_algo_config_id(func, theta);
        if ~isempty(algo_config_id)
            algorun_config_id = get_algorun_config_id(func.env.algo, algo_config_id, inst_id);
            if ~isempty(algorun_config_id)
                [runtime, censored, runlength, best_sol, solution_status] = get_runtime(algorun_config_id, seed, cutoff, dorunlength);
                if ~isempty(runtime)
                    return
                end
            end
        end
    end
end


%=== Quote a string: replace " by \" to run the command remotely.
function quoted = quote(string)
quoted = regexprep(string, '"', '\\"');

function str = random_string
s = rand('twister');
rn = rand;
rand('twister',s);
str = datestr(now, 'mmmm_dd_yyyy_HH_MM_SS_FFF');
str = strcat(str, '-', num2str(rn));