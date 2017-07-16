#=== Read inputs.
exec_path = ARGV[0]
command = ARGV[1...ARGV.length].join(" ")

#=== vsh comes with SecureCRT: C:\Program Files\SecureCRT\vsh.exe
cmd_connect = "vsh -accepthostkeys -auth publickey -noprompt -i z:\\.ssh\\id_rsa.pub hutter@arrow.cs.ubc.ca"

cmd_there = "cd #{exec_path}; #{command}"
cmd = cmd_connect + " " + cmd_there

#=== Need to read stdout. Got this tip from http://blade.nagaokaut.ac.jp/cgi-bin/scat.rb/ruby/ruby-talk/167193
result = nil
IO.popen(cmd + " 2>&1") {|file|
	result = file.gets(nil)
}
puts result