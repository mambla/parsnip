function [ ytest ] = go_tree()

% ClassificationTree upon the 8 SNPs before and after
% the one we are looking for
%prune to level 2
ytest = go_generic(tree_algorithm(8,0));

end