function [featureNdx, th, a , b, error] = selectBestRegressionStumpShare(goodi,cx, cz, cw);
% [th, a , b] = fitRegressionStump(x, z);
% z = a * (x>th) + b;
%
% where (a,b,th) are so that it minimizes the weighted error:
% error = sum(w * |z - (a*(x>th) + b)|^2) / sum(w)
%assumes cw is normalized to have a sum of 1

% atb, 2003
% torralba@ai.mit.edu

[Nfeatures, Nsamples] = size(cx{i}); 
numproblems = lenght(cx);

th = zeros(numproblems,Nfeatures);
a = zeros(numproblems,Nfeatures);
b = zeros(numproblems,Nfeatures);
error = zeros(numproblems,Nfeatures);

for m = 1:numproblems,
  tmpx = cx{m}(:,goodi);
  tmpz = cz{m}(goodi);
  tmpw = cw{m}(goodi);
  for n = 1:Nfeatures
    [th(m,n), a(m,n) , b(m,n), error(m,n)] = fitRegressionStump(tmpx(n,:), tmpz, tmpw);
  end
end

error = sum(error); %should this be square???

[error, featureNdx] = min(error);
th = th(:,featureNdx);
a = a(:,featureNdx);
b = b(:,featureNdx);
