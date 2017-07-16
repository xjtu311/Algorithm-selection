function rundata = dorun(func, theta_idxs, pi_idxs, captimes, seeds, iteration, rundata)
%=== Output which runs we're doing and do them.
assert(length(theta_idxs) == length(pi_idxs));
assert(length(theta_idxs) == length(seeds));
assert(length(theta_idxs) == length(captimes));

global runTimeForRunningAlgo;
runtimeBefore = runTimeForRunningAlgo;

for i=1:length(theta_idxs)
    bout(sprintf(strcat(['Iteration ', num2str(iteration), ': running config ', num2str(theta_idxs(i)) , ' on instance ', num2str(pi_idxs(i)), ' with seed ', num2str(seeds(i)), ' and captime ', num2str(captimes(i)),'\n'])));
end

tic
s=rand('twister');
[yNew, censoredNew, runtimeNew, runlengthNew, solvedNew, best_solNew] = func.funcHandle(theta_idxs, pi_idxs, seeds, captimes, func.tuningScenario);

if ~(func.cheap || func.matlab_fun)
    runtimeNew = max(runtimeNew,0.005);
end
rand('twister', s);
walltime = toc;

runtimeOfRun = max( 0.1, runTimeForRunningAlgo-runtimeBefore );
%=== Minimum counting 0.1 seconds per batch run to accout for overhead for
%=== calling the algorithm.
rundata.runTime(iteration+1) = rundata.runTime(iteration+1) + runtimeOfRun;
bout(sprintf(strcat(['Iteration ', num2str(iteration), ': doing ', num2str(length(theta_idxs)) , ' runs: took ', num2str(runTimeForRunningAlgo-runtimeBefore), 'CPU s, wall clock time: ', num2str(walltime), 's\n'])));
%alTunerTime = alTunerTime + (runTimeForRunningAlgo-runtimeBefore);

%=== When 100 nodes run jobs that each last milliseconds and request a DB connection things break (network goes down).
%=== Thus, sleep for a while so we don't create problems (of course, eventually, we'd want to do the whole interface via Matlab and
%=== just keep the one DB connection we have. (TODO)
minRuntimeToKeepNetworkTrafficDown = 0.1;
if (runTimeForRunningAlgo-runtimeBefore) < minRuntimeToKeepNetworkTrafficDown-1e-4 && func.db
    pause(minRuntimeToKeepNetworkTrafficDown - (runTimeForRunningAlgo-runtimeBefore))
end


for i=1:length(theta_idxs)
    %=== Check for an existing run with that theta_idx, pi_idx, and seed.
    idx = get_idx_for_theta_pi_seed(rundata, theta_idxs(i), pi_idxs(i), seeds(i));
    
    if ~isempty(idx)
        %=== Assert our current run has at least the same captime.
        assert(length(idx) == 1);
        orig_captime = rundata.used_captimes(idx);
        assert(captimes(i) >= orig_captime - 1e-6);
        
        %=== Delete the original entry, such that we don't have to 
        %=== worry about duplicates.
        all_but_idx = [1:idx-1, idx+1:length(rundata.used_theta_idxs)];
        rundata.used_theta_idxs = rundata.used_theta_idxs(all_but_idx);
        rundata.used_instance_idxs = rundata.used_instance_idxs(all_but_idx);
        rundata.usedSeeds = rundata.usedSeeds(all_but_idx);
        rundata.used_captimes = rundata.used_captimes(all_but_idx);
        rundata.y = rundata.y(all_but_idx);
        rundata.cens = rundata.cens(all_but_idx);
        rundata.runtimes = rundata.runtimes(all_but_idx);
        rundata.runlengths = rundata.runlengths(all_but_idx);
        rundata.solveds = rundata.solveds(all_but_idx);
        rundata.best_sols = rundata.best_sols(all_but_idx);
        rundata.iterations = rundata.iterations(all_but_idx);
        rundata.time_until_here = rundata.time_until_here(all_but_idx);
    end
    
    %=== Include new run in rundata.
    rundata.used_theta_idxs = [rundata.used_theta_idxs; theta_idxs(i)];
    rundata.used_instance_idxs = [rundata.used_instance_idxs; pi_idxs(i)];
    rundata.usedSeeds = [rundata.usedSeeds; seeds(i)];
    rundata.used_captimes = [rundata.used_captimes; captimes(i)];
    rundata.y = [rundata.y; yNew(i)];
    rundata.cens = [rundata.cens; censoredNew(i)];
    rundata.runtimes = [rundata.runtimes; runtimeNew(i)];
    rundata.runlengths = [rundata.runlengths; runlengthNew(i)];
    rundata.solveds = [rundata.solveds; solvedNew(i)];
    rundata.best_sols = [rundata.best_sols; best_solNew(i)];
    rundata.iterations = [rundata.iterations; iteration];
    rundata.time_until_here = [rundata.time_until_here; runTimeForRunningAlgo]; % same for all algo runs in the iteration
end

assert(all(~isnan(rundata.used_captimes)));
assert(all(~isnan(rundata.y)));
assert(all(~isnan(rundata.used_theta_idxs)));
assert(all(~isnan(rundata.usedSeeds)));