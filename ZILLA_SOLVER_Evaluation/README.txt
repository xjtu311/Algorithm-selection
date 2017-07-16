# This is a tool for building ZILLA based on cost sensitive DF.
# By doing so, you can evaluate the cost of omission for each solver

# Please check the path settings for all code since the intermediate results are solved locally.
# Please check script dir for further information

Key code:

For apply modeling and do solver subset selection (10 folds):

apply_zilla_models_SAT.m (better use compiled code)
apply_zilla_models_SAT_omiss.m (better use compiled code)


For building models for different presolver choices (10 folds):
build_zilla_models_SAT.m (better use compiled code)
build_zilla_models_SAT_omission.m (better use compiled code)


For obtain the final results for each fold:
testperformance_SAT_improve.m
testperformance_SAT_improve_omiss.m


For any question, please contanct Lin Xu at xulin730@cs.ubc.ca.


2-10-2012


