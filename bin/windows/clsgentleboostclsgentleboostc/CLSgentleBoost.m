function Model = CLSgentleBoost(X,y,sPARAMS);
%T requires sPARAMS.ux, sPARAMS.gamma, sPARAMS.startT

warning off
addpath boosting
warning on

if nargin<3
  sPARAMS.Nrounds = 100;
  sPARAMS.R = 0;
  sPARAMS.T = 0;
else
  if ~isfield(sPARAMS,'R');
    sPARAMS.R = 0;
  end
  if ~isfield(sPARAMS,'T');
    sPARAMS.T = 0;
  end
end

Nrounds = sPARAMS.Nrounds;
if length(Nrounds) == 1,
  Nrounds = [Nrounds Nrounds];
end

uniquey = unique(y);

warning off MATLAB:divideByZero
if ~sPARAMS.R & ~sPARAMS.T,
  Model.classifier = ...
      gentleBoostInsist(X,2*(y~=uniquey(1))'-1,Nrounds);
elseif sPARAMS.R
  if ~isfield(sPARAMS,'RECURSIVE'),
    sPARAMS.RECURSIVE = 0;
  end
  Model.classifier = ...
      gentleBoostR(X,2*(y~=uniquey(1))'-1,Nrounds(1),.1,sPARAMS.RECURSIVE);
elseif sPARAMS.T==1
  Model.classifier = ...
      gentleBoostT(X,2*(y~=uniquey(1))'-1,sPARAMS,Nrounds(1));
elseif sPARAMS.T==2
  [mclassifier,aa] = ...
      gentleBoostTbkg(X,2*(y~=uniquey(1))'-1,sPARAMS,Nrounds(1));  
  Model.classifier = mclassifier;
  Model.aa = aa;
else
  error('unknown R or T option');
end  
warning on MATLAB:divideByZero
fprintf('\n');

if ~iscell(Model.classifier)
  Model.Nbagging = 0;
  featuresused = zeros(Nrounds(1),1);
  for i =1:Nrounds, featuresused(i) = ...
      Model.classifier(i).featureNdx; end;
  Model.featuresused = featuresused;
else
  Nbagging = length(Model.classifier);
  Model.Nbagging = Nbagging;
  featuresused = zeros(Nrounds(1),Nbagging);
  for j = 1:Nbagging,for i =1:Nrounds, 
      featuresused(i,j) = Model.classifier{j}(i).featureNdx; 
  end; end
  Model.featuresused = featuresused;
end
 
Model.uniquey = uniquey;
Model.Nrounds = Nrounds;
Model.R = sPARAMS.R;
Model.T = sPARAMS.T;
