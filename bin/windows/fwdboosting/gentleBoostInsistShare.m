function classifier = gentleBoostInsistShare(cx, cy, Nrounds, beta)
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

numproblems = length(cx);
Nsamples = zeros(numproblems,1);
cw = cell(numproblems,1);

cclassifier = cell(numproblems,1);

for i = 1:numproblems,
  [Nfeatures, Nsamples(i)] = size(x); % Nsamples = Number of thresholds that we will consider
%  Fx = zeros(1, Nsamples);
  cw{i}  = ones(1, Nsamples); w = w/sum(w);
end

%for stopping criteria
chosen = zeros(Nfeatures,1);
m = 1;
pleasestop = 0;


while ~pleasestop %for m = 1:Nrounds
    fprintf('%d#',m);
    % weight triming
    ws = sort(w);
    cs = cumsum(ws); cs = cs/max(cs); 
    [foo, ndx] = min(abs(cs - beta));
    J = find(w >= ws(ndx));

    % weak regression 
    [featureNdx, vth, va , vb, verror] = selectBestRegressionStumpShare(J,cx,cy,cw);
    chosen(featureNdx) = chosen(featureNdx) + 1;
  
    % update parameters classifier
    for i = 1:numproblems,
      cclassifier{i}(m).featureNdx = featureNdx;
      cclassifier{i}(m).th = vth(i);
      cclassifier{i}(m).a  = va(i);
      cclassifier{i}(m).b  = vb(i);
      
      %%% Updating and computing classifier output
      fm = (va(i) * (cx{i}(featureNdx,:)>vth(i)) + vb(i));
      %Fx = Fx + fm;
      cw{i} = cw{i} .* exp(-cy{i}.*fm);
      cw{i} = cw{i} / sum(cw{i});
    end
    
    m = m + 1;
    pleasestop = ((sum(~~chosen)>Nrounds(1)) | (m > Nrounds(2)));
end
