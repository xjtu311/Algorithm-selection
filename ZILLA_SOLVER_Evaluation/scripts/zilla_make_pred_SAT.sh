#!/bin/tcsh
# setenv LD_LIBRARY_PATH /cs/local/generic/lib/pkg/matlab-7.8/bin/glnx86
setenv LD_LIBRARY_PATH /cs/local/generic/lib/pkg/matlab-7.10/runtime/glnx86:/cs/local/generic/lib/pkg/matlab-7.10/bin/glnx86
setenv DISPLAY fake_display
/ubc/cs/project/arrow/projects/FORLIN/ZILLA11/HAL_ZILLA/build_zilla_models_SAT $1 $2 $3
