function [ ytest ] = go_svm2()

% radial svm with g = 0.007 and c = 100 upon the 9 SNPs before and after
% the one we are looking for
% take the most 10 correlated permutated snps as features

ytest = go_generic(svm_algorithm2('-t 2 -g 0.007 -c 100 -q', 9, 10));

end