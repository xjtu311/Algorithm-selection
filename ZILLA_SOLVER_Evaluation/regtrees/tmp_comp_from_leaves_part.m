function [means, vars, tree_means] = tmp_comp_from_leaves_part(cell_of_leaves, numThetaIdx, numInstances, deterministic)
%ONLY FOR DEBUGGING the MEX CODE compute_from_leaves_part.c
M = length(cell_of_leaves);
tree_means = -ones(numThetaIdx,M);
tree_vars = -ones(numThetaIdx,M);

for m=1:M
    leaves = cell_of_leaves{m};

    tree_means(:,m) = 0;
    tree_vars(:,m) = 0;
    
    %=== Do the work for each leaf.
    for l_idx=1:length(leaves)
        leaf = leaves{l_idx};
        
        %=== Using mean & var from parametric fit.
        leaf_mean = leaf{5}(4);
        leaf_var  = leaf{5}(5);
        
        theta_idx_here = leaf{1};
        if deterministic
            error 'need to implement deterministic'
        else
            ratio_inst_idx = length(leaf{2})/numInstances;
            tree_means(theta_idx_here,m) = tree_means(theta_idx_here,m) + ratio_inst_idx * leaf_mean;
            tree_vars(theta_idx_here,m) = tree_vars(theta_idx_here,m) + ratio_inst_idx^2 * leaf_var;
        end
    end
end
means = mean(tree_means,2);
vars = mean(tree_vars + tree_means.^2,2) - means.^2;
assert(~any(isnan(vars)));
assert(all(vars>-1e-6));
