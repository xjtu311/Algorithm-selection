#! /usr/bin/X11/ruby
# you need change the path to make this runable
todofile=ARGV[0]
File.open(todofile){|lines|
  while oneline=lines.gets
    thisline=oneline.split(', ')
    dataset=thisline[0].chomp
    mymodel=thisline[1].chomp

shdir = "/ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/scripts/tmp"
         shfile="#{shdir}/qsub_#{dataset}_#{mymodel}_#{rand()}.sh"
         ff=File.new(shfile,'w')
         ff.puts "#! /bin/sh"
         ff.puts "#\$ -cwd"
         ff.puts "#\$ -o /ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/scripts/output  -e /ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/scripts/output"
         ff.puts "/ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/scripts/zilla_train_pre_SAT.sh #{dataset} $SGE_TASK_ID 1234 #{mymodel}"
         ff.close
                mycmd = "qsub -P eh -t 1-144 #{shfile}"
            puts mycmd
          system mycmd
end
}
