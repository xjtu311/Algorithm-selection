function runtime = get_algorun_result(func, x, inst_name, seed, cutoff, db)
%=== Returns the runtime of the algorithm on the specified instance with
%=== the specified seed and cutoff time.
%=== If using the database, get the appropriate IDs and check if the run is
%=== already stored there. 
%=== Otherwise, call a script via the command line that executes the run.
    if nargin < 5
        db = 0;
    end
    env = func.env;
    %=== Construct param_string.
    param_string = '';
    for i=1:length(func.cat)
        values = func.cat_all_values{i};
        param_string = strcat([param_string ' -'  func.cat_param_names{i} ' ' values{x(i)}]);
    end
    %=== Round continuous parameters to 4 digits -- otherwise DB problems !
    x(func.cont) = floor(x(func.cont).*10000 + 0.5)/10000;
    for i=1:length(func.cont)
        param_string = strcat([param_string ' -'  func.cont_param_names{i} ' ' num2str(x(i))]);
    end
            cmd = strcat(['ruby algoexec_wrapper.rb ' env.exec_path ' ' env.executable ' ' inst_name ' none ' num2str(cutoff) ' 2147483647 ' num2str(seed) ' ' param_string])

    %=== Get the result.
    if db
        %=== Check for entry in DB, if not there run the algo.
        runtime = ask_db(func, x, inst_name, cutoff, seed);
        if isempty(runtime)
            %=== Since the DB couldn't get us the value, start a run using the script runstarter.rb.
            run_cmd = strcat(['ruby runstarter.rb ' env.algo ' ' env.exec_path ' ' inst_name ' 0 ' num2str(cutoff) ' ' num2str(seed) ' ' param_string]);
            if isunix
                connect_cmd = strcat(['cd ' env.script_path '; ']); 
                cmd = strcat([connect_cmd ' ' run_cmd]);
                %[a,b] = unix(cmd);
                call(cmd); % faster version, but no outputs
            else
%                error('Dont want to run algos from DOS')
                connect_cmd = strcat(['ruby remote_executer.rb ' env.script_path]);
                cmd = strcat([connect_cmd ' ' run_cmd])
                [a,b] = dos(cmd);
            end
        end
        runtime = ask_db(func, x, inst_name, cutoff, seed);
        assert(~isempty(runtime), 'After executing the run, the runtime HAS to be in the DB.')
    else
        %=== Call command line script that returns the result once it's done.
        if isunix
            % faster version, but no outputs
            call(strcat(['ruby ~/arrowspace/altuning/scripts/algoexec_wrapper_for_linux.rb ' env.exec_path ' ' env.executable ' ' inst_name ' none ' num2str(cutoff) ' 2147483647 ' num2str(seed) ' ' param_string]));
%             [a,b] = unix(strcat(['ruby ~/arrowspace/altuning/scripts/algoexec_wrapper_for_linux.rb ' env.exec_path ' ' env.executable ' ' inst_name ' none ' num2str(cutoff) ' 2147483647 ' num2str(seed) ' ' param_string]))
%             b = b(64:end); % get rid of stupid extra output about terminal settings.
%             assert(~a);
%             runtime = str2num(b);
%             assert(runtime ~= -1);
        else
            cmd = strcat(['ruby algoexec_wrapper.rb ' env.exec_path ' ' env.executable ' ' inst_name ' none ' num2str(cutoff) ' 2147483647 ' num2str(seed) ' ' param_string])
            [a,b] = dos(cmd)
    %        [a,b] = unix('ruby ~/arrowspace/altuning/algoexec_wrapper_for_linux.rb ~/arrowspace/altuning/branin ruby branin_wrapper.rb dummy 0 5 100000000000 3 "-x 3 -y 3"')        
            assert(~a);
            runtime = str2num(b);
            assert(runtime ~= -1);
        end
    end
end
    
function runtime = ask_db(func, x, inst_name, cutoff, seed)
    runtime = [];
    inst_id = get_inst_id(inst_name);
    if ~isempty(inst_id)
        algo_config_id = get_algo_config_id(func, x);
        if ~isempty(algo_config_id)
            algorun_config_id = get_algorun_config_id(func.env.algo, algo_config_id, inst_id);
            if ~isempty(algorun_config_id)
                runtime = get_runtime(algorun_config_id, seed, cutoff);
            end
        end
    end
end