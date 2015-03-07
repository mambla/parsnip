function classifier = gentleBoostTree(x, y, Nrounds, beta)
% despite the name this is going to be ADABOOST...
%
% features x
% class: y = [-1,1]
%
% beta = weight triming. 
%
% Friedman, J. H., Hastie, T. and Tibshirani, R. 
% "Additive Logistic Regression: a Statistical View of Boosting." (Aug. 1998) 

% atb, 2003

if nargin < 4
    % default parameter for weight triming
    beta = .1;
end

[Nfeatures, Nsamples] = size(x); % Nsamples = Number of thresholds that we will consider

w  = ones(1, Nsamples); w = w/sum(w);

for m = 1:Nrounds
    fprintf('%dC',m);
    % weight triming
    ws = sort(w);
    cs = cumsum(ws); cs = cs/max(cs); 
    [foo, ndx] = min(abs(cs - beta));
    J = find(w >= ws(ndx)); 

    % weak regression 
    %[featureNdx, th, a , b, error] = select...Stump(x(:,J), y(J), w(J));
    syw.y = y(J);
    syw.w = w(J)';
    Tree=treefitweighted(x(:,J)',syw,'method','classification');
    % update parameters classifier
    classifier(m).Tree = Tree;
    
    % Updating and computing classifier output
    cry = Tree.classname(treeval(Tree,x'));
    for iii = 1:length(cry),
      ry(iii) = str2num(cry{iii});
    end
    corrects = (ry==y');
    et = sum((1-(corrects)).*w);
    alphat = .5*log((1-et)./et);
    classifier(m).alpha = alphat;
    
    w = w .* (exp(1-2*corrects));
    w = w / sum(w);
end

