function [ ytest ] = go_histogram()

% polynomial svm of degree=2 and c=1 is used upon descriptors containing
% histograms of all sequences of length 4 (sliding window interval is 1)
ytest = go_generic(histogram_algorithm('-t 1 -d 2 -c 1 -q',4,1));

% save result
mkdir('../results/histogram/');
save('../results/histogram/ytest.mat','ytest')

end