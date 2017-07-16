#! /usr/bin/X11/ruby
# you need change the path to make this work
dataset = ["SATRAND", "SATHAND", "SATINDU"]
shdir = "/ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/scripts/tmp"
for i in 0...dataset.length
         shfile="#{shdir}/qsub_#{dataset[i]}_#{rand()}.sh"
         ff=File.new(shfile,'w')
         ff.puts "#! /bin/sh"
         ff.puts "#\$ -cwd"
         ff.puts "#\$ -o /ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/scripts-11-16-2011/output  -e /ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/scripts-11-16-2011/output"
         ff.puts "/ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/scripts-11-16-2011/zilla_make_pred_SAT.sh #{dataset[i]} $SGE_TASK_ID 1234 "
         ff.close
                mycmd = "qsub -P eh -t 1-144 #{shfile}"
            puts mycmd
          system mycmd
end
