function [ys, censoreds, runtimes, runlengths, solveds, best_sols] = get_single_results(func, Theta_idx, instance_numbers, seeds, instance_filenames, censorLimits, kappa, dorunlength, name, numrun)

if ~func.db
    [ys, censoreds, runtimes, runlengths, solveds, best_sols] = get_single_results_nodb(func, Theta_idx, seeds, instance_filenames, censorLimits);
    return
end

%=== Returns the result (runtime/runlength) of the algorithm for the parameter configurations indexed by Theta_idx on the given instances with the given seeds.
%=== If using the database, get the appropriate IDs and check if the run is
%=== already stored there.
%=== Otherwise, call a script via the command line that executes the runs.

%TODO: clean up censoring for runlength

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

oldLeftTodo = 1:length(seeds);
assert(~isempty(oldLeftTodo));
numTry = 1;
index_to_run = [];
Theta_idx_to_run = [];
insts_to_run = {};
seeds_to_run = [];
censorLimits_to_run = [];
while ~isempty(oldLeftTodo)
    try % watch out for any exception from calling the target algorithm.
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
            end

            if (~func.db) || isempty(runtime)
                index_to_run(end+1) = i;
                Theta_idx_to_run(end+1,:) = theta_idx;
                insts_to_run{end+1} = instance_filename;
                seeds_to_run(end+1) = seed;
                if dorunlength
                    censorLimits_to_run(end+1) = 5;
                else
                    censorLimits_to_run(end+1) = censorLimit;
                end
            else
                runlengths(i) = runlength;
                runtimes(i) = runtime;
                censoreds(i) = censored;
                solveds(i) = solution_status;
                best_sols(i) = best_sol;

                leftTodo = setdiff(leftTodo, i);
            end

            if ~isempty(insts_to_run) && mod(length(insts_to_run), 10001)==0 && j < length(oldLeftTodo)
                % Start runs (don't wait) and empty queque of runs to do.
                [files_to_delete, leftTodo, solveds, censoreds, runtimes, runlengths, best_sols, Theta_idx_to_run, insts_to_run, seeds_to_run, censorLimits_to_run] = startRunsInBatch(func, Theta_idx_to_run, insts_to_run, seeds_to_run, censorLimits_to_run, index_to_run, name, leftTodo, solveds, censoreds, runtimes, runlengths, best_sols);
                delete_files(files_to_delete);
            end
        end

        if ~isempty(insts_to_run)
            [files_to_delete, leftTodo, solveds, censoreds, runtimes, runlengths, best_sols, Theta_idx_to_run, insts_to_run, seeds_to_run, censorLimits_to_run] = startRunsInBatch(func, Theta_idx_to_run, insts_to_run, seeds_to_run, censorLimits_to_run, index_to_run, name, leftTodo, solveds, censoreds, runtimes, runlengths, best_sols);
            delete_files(files_to_delete);
        end
        oldLeftTodo = leftTodo
        
    catch ME
        bout(strcat(['Num try: ', num2str(numTry)]));
        bout(ME.message);
        for i=1:length(ME.stack)
            bout(strcat(['In method ', ME.stack(i).file, ' at line ', num2str(ME.stack(i).line), '\n']));
        end
        
        numTry = numTry+1;
        if numTry >= 20
            func
%             Theta_idx
%             instance_filenames
%             seeds
%             censorLimits
%             name
%             oldLeftTodo
%             numTry
%             insts_to_run


            error('Still not successful after 20 tries.');
        end
        %Re-initialize these fields.
        index_to_run = [];
        Theta_idx_to_run = [];
        insts_to_run = {};
        seeds_to_run = [];
        censorLimits_to_run = [];
    end
end

global runTimeForRunningAlgo;
assert(all(runtimes >= -1e-5));
actual_runtimes = min([censorLimits'; runtimes'], [], 1); % to deal with inf in runs with wrong answer
runTimeForRunningAlgo = runTimeForRunningAlgo + sum(actual_runtimes);

switch func.singleRunObjective
    case 'runtime'
        ys = runtimes;
    case 'runlength'    
        ys = runlengths;
    case 'solqual'
        ys = best_sols;
        censoreds = zeros(size(ys,1),1);
    otherwise 
        error 'Still have to implement objective functions for single algorithm runs other than runtime, runlength, and solqual.'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% HELPER FUNCTIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delete_files(files_to_delete)
%=== Clean up files in directory RUNFILES - can't do earlier due to weird
%=== error (probably since Matlab is looking ahead, executing the delete
%=== too early?)
for i=1:length(files_to_delete)
    fprintf(strcat('Deleting file: ', files_to_delete{i}));
    delete(files_to_delete{i});
end

function [files_to_delete, leftTodo, solveds, censoreds, runtimes, runlengths, best_sols, Theta_idx_to_run, insts_to_run, seeds_to_run, censorLimits_to_run] = startRunsInBatch(func, Theta_idx_to_run, insts_to_run, seeds_to_run, censorLimits_to_run, index_to_run, name, leftTodo, solveds, censoreds, runtimes, runlengths, best_sols)
%=== This starts a batch of runs, and if we're not using the DB, parses the
%=== result and returns it. (If we use the DB, we read the result from
%=== there later.)
[files_to_delete, results] = startManyRuns(func, Theta_idx_to_run, insts_to_run, seeds_to_run, censorLimits_to_run, 1, name);
if ~func.db
    for k=1:length(index_to_run)
        i_prime = index_to_run(k);
        res = results(k,:);
        seed_out = res(end);    
        assert(seeds_to_run(k) == seed_out);

        solved = res(end-4);
        solveds(i_prime) = solved;
        if solved == 3 % wrong answer!
            warnstring = strcat(['Wrong answer for theta_idx ', num2str(Theta_idx_to_run(k)), ', instance ', num2str(insts_to_run(k)), ', and seed ', num2str(seeds_to_run(k)), '.\n']);
            warning(warnstring);
            bout(warnstring);

            censoreds(i_prime) = 1;
            runtimes(i_prime) = inf;
            runlengths(i_prime) = inf;
            best_sols(i_prime) = inf;
        else
            if solved == 0 % TIMEOUT
                censoreds(i_prime) = 1;
            else
                censoreds(i_prime) = 0;
            end
            runtimes(i_prime) = res(end-3);
            runlengths(i_prime) = res(end-2);
            best_sols(i_prime) = res(end-1);
        end
    end
    leftTodo = setdiff(leftTodo, index_to_run);
end

Theta_idx_to_run = [];
insts_to_run = {};
seeds_to_run = [];
censorLimits_to_run = [];



function [remote_files_to_delete, res] = startManyRuns(func, Theta_idx, instance_filenames, seeds, censorLimits, waitToFinish, name) 
remote_files_to_delete = [];
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
%deldir_post = strcat('DEL3RUNFILES/');
deldir_post = strcat('DEL4RUNFILES/');
if ~exist(strcat(func.outdir, '/', deldir_post), 'dir')
    mkdir(strcat(func.outdir, '/', deldir_post));
end

for openTry=1:100
    file_with_cmds = strcat(deldir_post,strrep(strrep(datestr(now), ' ', '-'), ':', '-'), '--', num2str(rand));
    remote_file_with_cmds = strcat(func.remote_outdir, file_with_cmds);
    file_with_cmds = strcat(func.outdir, file_with_cmds);
    
    cmds_fid = fopen(file_with_cmds,'w');
    if cmds_fid ~= -1 
        break
    end
    bout(strcat(['Try ', num2str(openTry), ' failed to open file: ', file_with_cmds, '\n']));
    if openTry == 100
        errstr = strcat(['Cannot open file to write callstrings to: ', file_with_cmds, '\n']);
        error(errstr);
    end
end
rand('state', s);

outfile = strcat(deldir_post, env.algo, '-', random_string(), '-out.txt');
remote_outfile = strcat(func.remote_outdir, outfile);
outfile = strcat(func.outdir, outfile);

%=== Run the commands in the file, if clusterize then on the cluster.
if func.db && clusterize
    if waitToFinish>0
        oncluster = 6;
    else
        oncluster = 8;
    end
else
    oncluster = 0;
end

global TestTheta;


if func.db
    fprintf(cmds_fid, strcat([env.algo, '\n', env.exec_path, '\n', env.params_filename, '\n', num2str(oncluster), '\n']));
%    fprintf(cmds_fid, strcat([env.algo, ' ', env.exec_path, ' ', env.params_filename, ' ', num2str(oncluster), '\n']));
else
    fprintf(cmds_fid, strcat([env.executable, '\n', env.exec_path, '\n', env.params_filename, '\n', num2str(-1), '\n']));
%    fprintf(cmds_fid, strcat(['\''', env.executable, '\''', ' ', env.exec_path, ' ', env.params_filename, ' ', num2str(oncluster), '\n']));
end
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
try
    fclose(cmds_fid); % throws error if we're deleting the file after.
catch ME
    bout(strcat(['Cannot close cmds_fid: ', num2str(cmds_fid), ' for file ', file_with_cmds]));
end
remote_files_to_delete = {file_with_cmds}; % remote_file_with_cmds

    
%=== If we're using the DB: just do the runs via the command line, read
%=== the results from the DB.
run_cmd = strcat(['ruby scripts/al_run_configs_in_file.rb ', remote_file_with_cmds]);
%run_cmd = strcat(['ruby scripts/al_run_configs_in_file_UBCSAT.rb ', remote_file_with_cmds]);
if ~func.db
    %=== If we don't use the database: also specify an outfile, where the
    %=== results will be stored and that we can parse later.
    %outfile = strcat(env.exec_path, '/', env.algo, '-', random_string(), '-out.txt');
    run_cmd = strcat([run_cmd, ' ', remote_outfile]);
end

if isunix
    connect_cmd = strcat(['cd ' env.script_path '; ']);
    cmd = strcat([connect_cmd ' ' run_cmd])
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

res = '';
if ~func.db
    res = csvread(outfile);
    remote_files_to_delete{end+1} = outfile;
end

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