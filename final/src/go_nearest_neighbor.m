function [ ytest ] = go_nearest_neighbor()

% Use 3-NN with data multiplied by guassian filter with sigma = 3 and
% euclidean distance metric
ytest = go_generic(nearest_neighbor_algorithm(3, 3, 'euclidean'));

% save result
mkdir('../results/nearest_neighbor/');
save('../results/nearest_neighbor/ytest.mat','ytest')

end