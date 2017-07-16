#=== Read inputs.
exec_path = ARGV[0]
#executable = ARGV[1]
#instance_relname = ARGV[2].to_f
#instance_specifics = ARGV[3].to_f
#cutoff_time = ARGV[4].to_f
#cutoff_length = ARGV[5].to_f
#seed = ARGV[6].to_f
#paramstring = ARGV[7...ARGV.length].join(" ")

command = ARGV[1...ARGV.length].join(" ")
#puts "Calling from path #{exec_path}: #{command}"

#=== vsh comes with SecureCRT: C:\Program Files\SecureCRT\vsh.exe
cmd_connect = "vsh -accepthostkeys -auth publickey -noprompt -i z:\\.ssh\\id_rsa.pub hutter@arrow01.cs.ubc.ca"

cmd_there = "cd #{exec_path}; #{command}"
cmd = cmd_connect + " " + cmd_there

#puts "Calling: #{cmd}"

#=== Need to read stdout. Got this tip from http://blade.nagaokaut.ac.jp/cgi-bin/scat.rb/ruby/ruby-talk/167193
runtime = -1;
censored = -1;
IO.popen(cmd + " 2>&1") {|file|
	while line = file.gets
#		puts line
# Result for ParamILS: #{solved}, #{runtime}, #{runlength}, #{best_sol}, #{seed}
		if line =~ /Result for ParamILS: (.*)/
			res = $1.chomp
			solved, runtime, runlength, best_sol, seed = res.split(",").map{|x| x.strip}
			if solved == "TIMEOUT"
				censored = 1;
			else
				censored = 0;
			end
		end
	end
}
puts runtime
#puts censored
#puts "#{runtime} #{censored}"