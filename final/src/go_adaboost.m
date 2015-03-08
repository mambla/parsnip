function [ ytest ] = go_tree()

% 10 iterations adaboost upon the 40 SNPs before and after
% the one we are looking for
% weak classifiers: snp in index i has value x?
ytest = go_generic(adaboost_algorithm(10,40));

end