function [cclass,lab,mmax,FM] = multiclass(X,y,Xtest,Nrounds);

if ~isempty(y);
  uniquey = unique(y);
  numc = length(uniquey);
  cclass = cell(numc,1);
  warning off MATLAB:divideByZero
  for i = 1:numc,
    if length(Nrounds) == 1,
      cclass{i} = gentleBoost(X, 2*(y==uniquey(i))'-1, Nrounds);
    elseif length(Nrounds)==2
      cclass{i} = gentleBoostInsist(X, 2*(y==uniquey(i))'-1, ...
	  Nrounds);
    else 
      cclass{i} = gentleBoostR(X, 2*(y==uniquey(i))'-1, ...
	  Nrounds(1));
    end
  end
  warning on MATLAB:divideByZero
else
  cclass = X;
  numc = length(cclass);
  uniquey = Nrounds;
end

FM = zeros(size(Xtest,2),numc);
if nargin>2
  for i = 1:numc,
    [Cx,Fx] = strongGentleClassifier(Xtest,cclass{i});
    FM(:,i) = Fx';
  end
  [mmax,lab] = max(FM,[],2);
end

lab = uniquey(lab);
