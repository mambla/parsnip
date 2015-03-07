function classifier = gentleBoostR(x, y, Nrounds, beta, RECURSIVE)
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

if nargin<5
  RECURSIVE = 0;
end

[Nfeatures, Nsamples] = size(x); % Nsamples = Number of thresholds that we will consider

if RECURSIVE>0
  global TESTX;
  global TESTY;
  allfeaturesused = [];
  for i = 1:RECURSIVE,
    fprintf('(%d)',i);
    cclassifier{i} = gentleBoostR(x, y, Nrounds, beta, 0);
    for j = 1:Nrounds,
      allfeaturesused = [allfeaturesused;cclassifier{i}(j).featureNdx];
    end
  end
  if 0 
    was10 = 50;
    xx = zeros(Nfeatures,Nsamples*was10);
    yy = zeros(Nsamples*was10,1);
    nallfeaturesused = size(allfeaturesused,1);
    
    for i = 1:Nsamples*was10,
      featureNdx = floor(rand*allfeaturesused)+1;
      [tmpsort,randind] = sort(rand(Nsamples,1));
      xx(:,i) = x(:,randind(1));
      xx(featureNdx,i) = x(featureNdx,randind(2));
      yy(i) = y(randind(1));
    end
    for i = 1:RECURSIVE,
      ry = strongGentleClassifier(xx,cclassifier{i});
      trainingerror(i) = mean(ry'~=yy);
      if ~isempty(TESTX),
	rTESTY = strongGentleClassifier(TESTX,cclassifier{i});
	testingerror(i) = mean(rTESTY'~=TESTY);
      end
    end
    if ~isempty(TESTX),
      figure;
      plot(testingerror,trainingerror,'.');
      drawnow;
    end
    
    [minn,mini] = min(trainingerror);
    classifier = cclassifier{i};
  else
    classifier = cclassifier;
  end
  return;
end %end of RECURSIVE part

if nargin < 4
    % default parameter for weight triming
    beta = .1;
end

Fx = zeros(1, Nsamples);
w  = ones(1, Nsamples); w = w/sum(w);

for m = 1:Nrounds
  if ~mod(m,10),
      fprintf('%d',m);
    else
      fprintf('R');
    end
    % weight triming
    ws = sort(w);
    cs = cumsum(ws); cs = cs/max(cs); 
    [foo, ndx] = min(abs(cs - beta));
    J = find(w >= ws(ndx)); 

    % weak regression 
    [featureNdx, th, a , b, error] = selectBestRegressionStump(x(:,J), y(J), w(J));
    
    % update parameters classifier
    classifier(m).featureNdx = featureNdx;
    classifier(m).th = th;
    classifier(m).a  = a;
    classifier(m).b  = b;
    
    % Updating and computing classifier output
    fm = (a * (x(featureNdx,:)>th) + b);
    [tmpsort,randind] = sort(rand(Nsamples,1));
    x(:,end+1) = x(:,randind(1)); %also try from a different label,
                                  %or mean
    x(featureNdx,end) = x(featureNdx,randind(2));
    [tmpy,Fx(end+1)] = ...
	strongGentleClassifier(x(:,end),classifier(1:m));
    [tmpy,tmpFx] = ...
	strongGentleClassifier(x(:,end),classifier(1:(m-1)));
    fm(end+1) = Fx(end) - tmpFx;
    Fx(1:(end-1)) = Fx(1:(end-1)) + fm(1:(end-1));
    w(end+1) = w(randind(1));
    y(end+1) = y(randind(1));
    w = w .* exp(-y.*fm); %this is only approximated
    w = w / sum(w);
end

