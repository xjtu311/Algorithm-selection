function myassert(condition, errorString)
if nargin < 2
    errorString = 'Assertion failed.';
end
if ~condition
    sprintf(strcat('ERROR: ', errorString));
    error(errorString)
end
