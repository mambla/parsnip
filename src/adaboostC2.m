
function [ y_classifier_result, r ] = adaboostC2( X, Model )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% this run the classifier on the test samples.
% return as well r- how mutch the classifier is sure about it's answer

i_feature = Model.i;
thresh = Model.t;
p = Model.p;
alpha = Model.alpha;

% num of samples
n = size(X,2);

r = zeros(n,1);
y_classifier_result = zeros(n,1);

%for each sample
for k=1:n
    r(k) = sum(alpha .* ((2*(X(i_feature,k) == thresh)-1) .* p));
    y_classifier_result(k) = sign(r(k));
end

end
