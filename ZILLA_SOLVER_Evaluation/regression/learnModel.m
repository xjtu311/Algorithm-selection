function model = learnModel(X, y, cens, cat, catDomains, isclean, al_opts, names, varargin)
if nargin < 7
    names = {};
end
%=== Wrapper around initializing a model and then preparing it for prediction.
model = initModel(X, y, cens, cat, catDomains, isclean, al_opts, names);
if isfield(al_opts, 'params_to_use')
    model.params = al_opts.params_to_use;
else
    model = optimizeModel(model, varargin{:});
end
model = prepareModel(model);