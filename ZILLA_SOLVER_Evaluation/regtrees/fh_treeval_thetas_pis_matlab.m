function cell_of_leaves = fh_treeval_thetas_pis_matlab(Tree,Theta,X)
%   FH_TREEVAL_THETAS_PIS 
%   Propagate the parameter configruations in Theta's rows and the 
%   instances with the features in X's row down the tree.
%   Return a cell array with entries [theta_idx, pi_idx, value, nodenumber]
%   for each leaf that is compatible with a subset of both Theta and X.
%
%   NaN values in Theta or X are forbidden.

if ~isstruct(Tree) | ~isfield(Tree,'method')
   error('stats:treeval:BadTree',...
         'The first argument must be a decision tree.');
end
[N,nt] = size(Theta);
[n,nx] = size(X);

if nt+nx ~= Tree.npred
   error('stats:treeval:BadInput',...
         'The X and Theta matrices must have %d columns together.',Tree.npred);
end

cell_of_leaves = fwd(Tree,Theta,X,1:N,1:n,1);

%------------------------------------------------
function cell_of_leaves = fwd(regTree,Theta,X,Theta_rows,x_rows,thisnode)
%DOAPPLY Apply classification rule to specified rows starting at a node.
%   This is a recursive function.  Starts at top node, then recurses over
%   child nodes.  THISNODE is the current node at each step.

splitvar      = regTree.var(thisnode);
cutoff        = regTree.cut(thisnode);
assignedclass = regTree.class(thisnode);
kids          = regTree.children(thisnode,:);
catsplit      = regTree.catsplit;

% Terminal case
if splitvar==0
   id = regTree.class(thisnode);
   cell_of_leaves = {{[Theta_rows], [x_rows], id, thisnode}};
   return;
end

%%%% Now deal with non-terminal nodes %%%%
cell_of_leaves = {};

% Determine whether splitting on a parameter or an instance feature.
if abs(splitvar) <= size(Theta,2)
    split_on_param = 1;
    x = Theta(Theta_rows,abs(splitvar));   
else 
    split_on_param = 0;
    x = X(x_rows,abs(splitvar)-size(Theta,2));   
end

% Determine if this point goes left, goes right, or stays here
if splitvar>0                % continuous variable
  isleft = (x < cutoff);
  isright = ~isleft;
  ismissing = isnan(x);
else                         % categorical variable
  isleft = ismember(x,catsplit{cutoff,1});
  isright = ismember(x,catsplit{cutoff,2});
  ismissing = ~(isleft | isright);
end

subrows = find(isleft & ~ismissing);  % left child node
if ~isempty(subrows)
    if split_on_param
        cell_of_leaves = fwd(regTree,Theta,X,Theta_rows(subrows),x_rows,kids(1));
    else
        cell_of_leaves = fwd(regTree,Theta,X,Theta_rows,x_rows(subrows),kids(1));
    end
end

subrows = find(isright & ~ismissing); % right child node
if ~isempty(subrows)
    if split_on_param
        cell_right = fwd(regTree,Theta,X,Theta_rows(subrows),x_rows,kids(2));
    else
        cell_right = fwd(regTree,Theta,X,Theta_rows,x_rows(subrows),kids(2));
    end
    cell_of_leaves = [cell_of_leaves cell_right]; % concatenate into one cell array
end

subrows = find(ismissing);            % missing, treat as leaf.
if ~isempty(subrows)
    if split_on_param
        cell_of_quasi_leaves = {{Theta_rows(subrows), x_rows, regTree.class(thisnode), thisnode}};
    else
        error 'No empty features (continuous) allowed - a value is either bigger or smaller/equal another one.'
%        cell_of_quasi_leaves = {{Theta_rows, x_rows(subrows), regTree.class(thisnode)}};
    end
    cell_of_leaves = [cell_of_leaves cell_of_quasi_leaves]; % concatenate into one cell array
end