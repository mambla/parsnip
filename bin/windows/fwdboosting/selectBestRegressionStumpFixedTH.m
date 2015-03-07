function [featureNdx, th, a , b, error] = selectBestRegressionStumpFixedTH(x, z, w);
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
%    [th(n), a(n) , b(n), error(n)] = fitRegressionStump(x(n,:), z,
%    w);
     %[a(n),b(n)] = fitregression(x(n,:)>=0,z,w);
     b(n) = 0;
     %a*x = z
     a(n) = (w.*(x(n,:)))'\((w.*z)');
     errv = (a(n).*(x(n,:))+b(n)-z).^2;
     error(n) = sum(w.*errv);
end

[error, featureNdx] = min(error);
th = 0;
a = a(featureNdx);
b = b(featureNdx);
