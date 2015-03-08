function [ ytest ] = go_const()

% algorithm that guesses 1 for all SNPs
ytest = go_generic(template_algorithm(1));

% save result
mkdir('../results/const/');
save('../results/const/ytest.mat','ytest')

end