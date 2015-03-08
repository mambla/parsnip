function [ ytest ] = go_svm2()

% radial svm with g = 0.01 and c = 100 upon the 30 SNPs before and after
% the one we are looking for
% take the most 10 correlated permutated snps as features

ytest = go_generic(svm_algorithm2('-t 2 -g 0.01 -c 100 -q', 30, 10));

end