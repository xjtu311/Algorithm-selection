Four Steps for runing ZILLA solver evaluation:
Note: you may need change the path to make it work

1. ruby  subjob_zilla_pre_SAT.rb
 make prediction and write into files under SATMODELS

2. ruby subjob_train_perform_SAT.rb todofile
 do solver subset selection and obtain training performance. Write results into files

3.  testperformance_SAT_improve.m
 puts all folds together and generate the final results files

4. finalSAT.m
  analysis


