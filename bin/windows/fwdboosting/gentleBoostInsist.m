function classifier = gentleBoostInsist(x, y, Nrounds, beta)
% gentleBoost
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
Fx = zeros(1, Nsamples);

%EDITED: DELETED - original w was uniformly
%w  = ones(1, Nsamples); w = w/sum(w);

%EDITED: ADDED - my new initial w
w  = zeros(1, Nsamples);
pos_w = 0.5 * 1/sum(y==1);
neg_w = 0.5 * 1/sum(y==-1);
w( y==1 )= pos_w; 
w( y==-1 )= neg_w; 
%end EDITION

chosen = zeros(Nfeatures,1);
m = 1;
pleasestop = 0;
while ~pleasestop %for m = 1:Nrounds
    if ~mod(m,10),
      fprintf('%d',m);
    else
      fprintf('#');
    end
    % weight triming
    ws = sort(w);
    cs = cumsum(ws); cs = cs/max(cs); 
    [foo, ndx] = min(abs(cs - beta));
    J = find(w >= ws(ndx)); 

    % weak regression 
    [featureNdx, th, a , b, error] = selectBestRegressionStump(x(:,J), y(J), w(J));
    chosen(featureNdx) = chosen(featureNdx) + 1;
  
    % update parameters classifier
    classifier(m).featureNdx = featureNdx;
    classifier(m).th = th;
    classifier(m).a  = a;
    classifier(m).b  = b;
    
    % Updating and computing classifier output
    fm = (a * (x(featureNdx,:)>th) + b);
    Fx = Fx + fm;
    w = w .* exp(-y.*fm);
    w = w / sum(w);
    m = m + 1;
    pleasestop = ((sum(~~chosen)>Nrounds(1)) | (m > Nrounds(2)));
end
