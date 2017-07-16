#!/bin/bash

scenario=$1
name=$2
aclib=$3

# Get parent directory of file
parentDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd $aclib

$parentDir/../../bin/hydra --num-iterations 8 --num-smac-runs 4 --num-configs-per-iter 2 --rungroup $name --num-run 1 --smacOptions $scenario --smac-execution-options $parentDir/smac-execution-options-local.txt --zilla-options $parentDir/zilla-options.txt --portfolio-evaluation ZILLA --experiment-dir $aclib