function Model = CLSadaBoost(X,y,inparams);
 
if nargin<3
  inparams.Nrounds = 100;
end

Model.inparams = inparams;
%Model.X = X;

if (isfield(inparams,'KO'))
  KO = inparams.KO;
else
  KO = 0;
end

if (isfield(inparams,'NoRepeats'))
  NoRepeats = logical(inparams.NoRepeats);
else
  NoRepeats = false;
end

if (isfield(inparams,'AdaBoostVC_d'))
  d = inparams.AdaBoostVC_d;
else
  d = 0;
end

if (not(isfield(inparams,'feature_limit')))
  inparams.feature_limit = +inf;
end

[nfeatures,npoints] = size(X);

if (NoRepeats)
  usable_features = 1:nfeatures;
end

alpha = zeros(inparams.Nrounds,1);
%W = zeros(npoints,inparams.Nrounds);
w = ones(npoints,1)/npoints;

chosen = zeros(nfeatures,1);

%cModels = cell(inparams.Nrounds,1);
for i = 1:inparams.Nrounds
  %fprintf(1,'@');
  w = w ./ sum(w);
  %W(:,i) = w;
  
  if (NoRepeats)
    [abbrev_featureNdx, th, error,a,b] = selectBestStump(X(usable_features,:), y, w);
    featureNdx = usable_features(abbrev_featureNdx);
    usable_features(abbrev_featureNdx) = [ ];
  else
    [featureNdx, th, error,a,b] = selectBestStump(X, y, w);  
  end
  
  chosen(featureNdx) = chosen(featureNdx) + 1;

  % update parameters classifier
  classifier.featureNdx = featureNdx;
  classifier.th = th;
  classifier.a = a;
  classifier.b = b;
  cModels{i} = classifier;
  
  if (KO)
    % Perform Feature Knockout:
    randind1 = ceil(rand * npoints);
    randind2 = ceil(rand * npoints);
    X(:,end+1) = X(:,randind1);
    X(featureNdx,end) = X(featureNdx,randind2);
    y(  end+1) = y(  randind1);
    w(  end+1) = w(  randind1);
  end

  [ry ,weight2] = adaBoostClassifier(X,classifier);
  
  e = sum(w.*(ry~=y));
  e = max(e,1/(100*npoints));
  
  if (d > 0)
    % AdaBoost-VC: Long, Vega, "Boosting and Microarray Data"
    m = npoints;
    e = e + d/m*(log(m)+(1 + e * m / d) .^ 0.5); % e = e_emp
  end
  
  Model.ve(i) = e;
  alpha(i) = 0.5*log((1-e)/e);
  w = w.*exp(-alpha(i)*(y.*ry));
  
  if (KO)
    fprintf('O');
  elseif (NoRepeats)
    fprintf('1');
  elseif (d>0)
    fprintf(sprintf('VC%d ', d));
  else
    fprintf('A');
  end
  
  if (nnz(chosen) >= inparams.feature_limit)
    fprintf('\nFeature limit hit at %d features!\n', inparams.feature_limit);
    break;
  end
  
  if (NoRepeats & isempty(usable_features))
    fprintf('\nNo features left to use!\n');
    break;
  end
   
end

Model.cModels = cModels;
Model.alpha = alpha;

Model.featuresused = find(chosen>0);

fprintf('\n');

