function [ ytest ] = go_svm()

% radial svm with g = 0.007 and c = 100 upon the 9 SNPs before and after
% the one we are looking for
ytest = go_generic(svm_algorithm('-t 2 -g 0.007 -c 100 -q', 9));

end