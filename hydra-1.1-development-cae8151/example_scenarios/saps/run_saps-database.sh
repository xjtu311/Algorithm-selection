#!/bin/bash
#PBS -S bin/bash

# some arguments don't have a corresponding value to go with it such
# as in the --default example).
while [[ $# > 0 ]]
do
key="$1"

case $key in
    -p|--pool)
    POOL="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

../../bin/hydra --num-iterations 4 --num-smac-runs 4 --num-configs-per-iter 2 --rungroup Hydra_Saps --num-run 1 --smacOptions ./saps-scenario.txt --smac-execution-options ./smac-execution-options-database.txt --zilla-options ./zilla-scenario.txt --portfolio-evaluation ZILLA --mysqldbtae-pool $POOL