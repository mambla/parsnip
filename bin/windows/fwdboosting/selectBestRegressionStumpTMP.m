function [featureNdx, th, a , b, error] = selectBestRegressionStumpTMP(x, z, w);
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
warning off
tic
for n = 1:Nfeatures
    [th(n), a(n) , b(n), error(n)] = fitRegressionStump(x(n,:), z, w);
end
toc
tic  
[thv,av,bv,errorv] = fitRegressionStumpMATRIX(x,z,w);
toc
warning on
norm(th'-thv)
norm(a'-av)
norm(b'-bv)
norm(error'-errorv)
[error, featureNdx] = min(error);
th = th(featureNdx);
a = a(featureNdx);
b = b(featureNdx);
