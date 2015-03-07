function [vth, va , vb, error,avec,bvec] = fitRegressionStumpMATRIX(x, z, w);
% [th, a , b] = fitRegressionStump(x, z, w);
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

[Nfeatures, Nsamples] = size(x); % Nsamples = Number of thresholds that we will consider
%%w = w/sum(w); % just in case...

[x, j] = sort(x,2); % this now becomes the thresholds
th = x;
z = z(j); w = w(j);

Szw = cumsum(z.*w,2); Ezw = Szw(:,end);
Ezwmat = Ezw*ones(1,Nsamples);
Sw  = cumsum(w,2);
%Sw = [0 Sw];
% This is 'a' and 'b' for all posible thresholds:
b = Szw ./ Sw;
a = (Ezwmat - Szw) ./ (1-Sw) - b;

% Now, let's look at the error so that we pick the minimum:
% the error at each threshold is:
% for i=1:Nsamples
%     error2(i) = sum(w.*(z - ( a(i)*(x>th(i)) + b(i)) ).^2);
%alt  error2(i) = sum(w.*(z ~= (x>th(i));
% end
% but with vectorized code it is much faster but also more obscure code:

%%error = sum(w.*z.^2) - ...
%%        2*a.*(Ezw-Szw) - 2*b*Ezw + (a.^2 +2*a.*b) .* (1-Sw) + b.^2;
error = sum(w.*z.^2,2)*ones(1,Nsamples) - 2*a.*(Ezwmat-Szw) - ...
    2*b.*(Ezwmat) + (a.^2 +2*a.*b) .* (1-Sw) + b.^2;
%alt: error = sum(w.*z.^2,2)*ones(1,Nsamples) - 2*(Ezwmat-Szw) + (1-Sw);

%[error; error2]
% Output parameters. Search for best threshold (th):
lengthe = size(error,2);
[error, k] = min(error,[],2);

if 0 
  th(:,end+1) = th(:,end);
  i = [(1:Nfeatures)';(1:Nfeatures)'];
  j = [k;(k+1)];
  s = ones(Nfeatures*2,1)*0.5;
  tmpmat = sparse(i,j,s);
  th = diag(th*tmpmat);
  i = (1:Nfeatures)';
  j = k;
  s = ones(Nfeatures,1);
  tmpmat = sparse(i,j,s);
  a = diag(a*tmpmat);
  b = diag(b*tmpmat);
elseif 0
  th(:,end+1) = th(:,end);
  va = zeros(Nfeatures,1);
  vb = zeros(Nfeatures,1);
  vth = zeros(Nfeatures,1);
  for i = 1:Nfeatures,
    vth(i) = (th(i,k(i)) + th(i,k(i)+1))/2;
    va(i) = a(i,k(i));
    vb(i) = b(i,k(i));
  end
else
  ind = sub2ind([Nfeatures,Nsamples],(1:Nfeatures)',k);
  va = a(ind);
  vb = b(ind);
  ind2 = sub2ind([Nfeatures,Nsamples],(1:Nfeatures)',min(k+1,Nsamples));
  vth = (th(ind) + th(ind2))/2;
end
%if k == lengthe  % can probably use length(th)
%    th = th(k);
%else
%    th = (th(k) + th(k+1))/2;
%end

%a = a(k);
%b = b(k);
