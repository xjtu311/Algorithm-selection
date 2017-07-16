function stopflag = have_to_stop(alTunerTime, runTime, func, used_instance_idxs, al_opts, iteration)
stopflag = ((alTunerTime + sum(runTime) > func.tuningTime) || (size(used_instance_idxs,1) >= al_opts.totalNumRunLimit)) || (iteration >= al_opts.numIterations);
%stopflag = ((sum(runTime) > func.tuningTime) || (size(used_instance_idxs,1) >= al_opts.totalNumRunLimit)) || (iteration >= al_opts.numIterations);
