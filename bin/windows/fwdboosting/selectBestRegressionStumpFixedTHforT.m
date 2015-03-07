function [featureNdx, th, a , b, error] = selectBestRegressionStumpFixedTHforT(x, z, w, select);
% [th, a , b] = fitRegressionStump(x, z);
% z = a * (x>th) + b;
%
% where (a,b,th) are so that it minimizes the weighted error:
% error = sum(w * |z - (a*(x>th) + b)|^2) / sum(w)

% atb, 2003
% torralba@ai.mit.edu

[Nfeatures, Nsamples] = size(x); % Nsamples = Number of thresholds that we will consider
w = w/sum(w); % just in case...

th = zeros(1,Nfeatures);
a = zeros(1,Nfeatures);
b = zeros(1,Nfeatures);
error = zeros(1,Nfeatures);

for n = 1:Nfeatures
  [tmpth, tmpa , tmpb, tmperror, av,bv,errorv] = fitRegressionStump(x(n,:), z, w);
  a(n) = av(select(n));
  b(n) = bv(select(n));
  error(n) = errorv(select(n));
end

[error, featureNdx] = min(error);
th = 0;
a = a(featureNdx);
b = b(featureNdx);
