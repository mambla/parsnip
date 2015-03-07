function [th, a , b, error,avec,bvec] = fitRegressionStump2(x, z, w);
% [th, a , b] = fitRegressionStump(x, z);
% The regression has the form:
% z = a * (x>th) + b;
%
% where (a,b,th) are so that it minimizes the weighted error:
% error = sum(w * |z - (a*(x>th) + b)|^2) 
%
% x,z and w are vectors of the same length
% x, and z are real values.
% w is a weight of positive values. There is no asumption that it sums to
% one.
% atb, 2003

[Nfeatures, Nsamples] = size(x) % Nsamples = Number of thresholds that we will consider
%Nsamples = length(x); % Nsamples = Number of thresholds that we will consider

w = w/sum(w); % just in case...

[x, j] = sort(x); % this now becomes the thresholds
th = x;
z = z(j); w = w(j);

%Sw = 1;
Swx = cumsum(w);
Swz = sum(w.*z);
Swxz = cumsum(z.*w);
%Swx2 = Swx
Swz2 = sum((z.^2).*w)
%Ezw = Szw(end);

num =  Swxz - Swx * Swz;
del1 = Swx - Swx.^2;
%del2 = Swz2 - Swz^2;

% This is 'a' and 'b' for all posible thresholds:
a = num ./ del1
b = (Swx*Swz - Swx.*Swxz)./del1;

% Now, let's look at the error so that we pick the minimum:
% the error at each threshold is:
%Nsamples
 for i=1:Nsamples
   [tmpa,tmpb] = fitregression(x>th(i),z,w);
   a(i) = tmpa;
   b(i) = tmpb;
     error2(i) = sum(w.*(z - ( a(i)*(x>th(i)) + b(i)) ).^2);
 end
% but with vectorized code it is much faster but also more obscure code:
% error = sum(w.*z.^2) - 2*a.*(Ezw-Szw) - 2*b*Ezw + (a.^2 +2*a.*b) .* (1-Sw) + b.^2;
%[error; error2]
% Output parameters. Search for best threshold (th):
error = error2
lengthe = length(error);
[error, k] = min(error);

if k == lengthe  % can probably use length(th)
    th = th(k);
else
    th = (th(k) + th(k+1))/2;
end
avec = a;
bvec = b;
a = a(k);
b = b(k);
