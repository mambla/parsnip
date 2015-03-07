function [Cx, Fx] = strongGentleClassifierTree(x, classifier)
% [Cx, Fx] = strongLogitClassifier(x, classifier)
%
% Cx is the predicted class 
% Fx is the output of the additive model
% Cx = sign(Fx)
%
% In general, Fx is more useful than Cx.
%
% The weak classifiers are stumps

% Friedman, J. H., Hastie, T. and Tibshirani, R. 
% "Additive Logistic Regression: a Statistical View of Boosting." (Aug. 1998) 

% atb, 2003
% torralba@ai.mit.edu

Nstages = length(classifier);
[Nfeatures, Nsamples] = size(x); % Nsamples = Number of thresholds that we will consider
numclasses = length(classifier(1).Tree.classname);
Fx = zeros(numclasses, Nsamples);
for m = 1:Nstages
  sfit = treeval(classifier(m).Tree,x');      % find assigned class numbers
  sfit = classifier(m).Tree.classname(sfit);    % get class names
  ry = zeros(size(sfit));
  for iii = 1:length(sfit),
    ry(iii) = str2num(sfit{iii});
  end
  for i = 1:numclasses,
    Fx(i,:) = Fx(i,:) + (classifier(m).alpha .* (ry==i))';
  end
end

[maxx,Cx] = max(Fx);
