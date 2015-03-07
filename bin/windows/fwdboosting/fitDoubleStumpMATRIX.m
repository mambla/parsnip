function [vth, va , vb, error,avec,bvec] = fitDoubleStumpMATRIX(x1,x2, z, w);
% This is not really necessary just pass max x1,x2 to regular code
% and also -min(x1,x2)
% [th, a , b] = fitRegressionStump(x, z, w);
% The regression has the form:
% z = a * (x1>th & x2>th) + b;
%
% where (a,b,th) are so that it minimizes the weighted error:
% error = sum(w * |z - (a*(x>th) + b)|^2) 
%
% x1,x2, are matrices. z and w are vectors of the same length
% x, and z are real values.
% w is a weight of positive values. There is no asumption that it sums to
% one.
% atb, 2003

[Nfeatures, Nsamples] = size(x1); 
% Nsamples = Number of thresholds that we will consider

x = max(x1,x2);

[x, j] = sort(x,2); % this now becomes the thresholds
th = x;
z = z(j); w = w(j);

Szw = cumsum(z.*w,2); 
Ezw = Szw(:,end);
Ezwmat = Ezw*ones(1,Nsamples);
Sw  = cumsum(w,2);
%Sw = [0 Sw];
% This is 'a' and 'b' for all posible thresholds:
b = Szw ./ Sw;
a = (Ezwmat - Szw) ./ (1-Sw) - b;

% Now, let's look at the error so that we pick the minimum:
% the error at each threshold is:

error = sum(w.*z.^2,2)*ones(1,Nsamples) - 2*a.*(Ezwmat-Szw) - ...
    2*b.*(Ezwmat) + (a.^2 +2*a.*b) .* (1-Sw) + b.^2;

% Output parameters. Search for best threshold (th):
lengthe = size(error,2);
[error, k] = min(error,[],2);

ind = sub2ind([Nfeatures,Nsamples],(1:Nfeatures)',k);
va = a(ind);
vb = b(ind);
ind2 = sub2ind([Nfeatures,Nsamples],(1:Nfeatures)',min(k+1,Nsamples));
vth = (th(ind) + th(ind2))/2;

