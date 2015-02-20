function [ Model ] = TrainAdaboost2( X, y, target, nRounds)
%X's column is sample. c's row is feature.
% y labels vector (0/1/2)
% target=0/1/2. we would build classifier that detect this value or not
% this value
%nRounds - num of adaboost rounds

sPARAMS.T = nRounds;

switch target
    case 0
        y = changem(y,[1 -1 -1], [0 1 2]);
    case 1
        y = changem(y,[-1 1 -1], [0 1 2]);
    case 2
        y = changem(y,[-1 -1 1], [0 1 2]);
end

Model = myAdaboost(X,y,sPARAMS);

end


function [ Model ] = myAdaboost( X, y, sParams)

if isfield(sParams,'T')
    T = sParams.T;
else
    T = 10;
end

i = ones(T,1);
theta = zeros(T,1);
p = zeros(T,1);
alpha = zeros(T,1);

% num of samples
n = size(X,2);
%initialize the weight of each sample
w = ones(n,1) / n;

% do T itererations of finding best classifier
for t = 1:T
   w = w./sum(w);
   [i(t),theta(t),p(t),e,y_guess] = selectwc(X,y,w);
   alpha(t) = 0.5 * log((1-e)/e);
   
   if (isinf(alpha(t)) || isnan(alpha(t)))
       alpha(t) = 0;
        break
    end
   
   %update w (weights)
   w = w .* exp(-1 * y .* alpha(t) .* y_guess);   
end

Model.i = i;
Model.t = theta;
Model.p = p;
Model.alpha = alpha;

end

function [ i, theta, p, e, guess_y ] = selectwc( X, y, w )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% d = num of attributes
% n = num of samples in the training data
[d, n] = size(X);
e = 2; % max Error possible

% for each attribute 
for index = 1:d    
    %for each decision stumb (options)
    for t = 0:2
        for pTmp=-1:2:1
          thresh = t;
          [tmp_e, tmp_y_guess] = run_weak_classifier(X,y,w,index,thresh,pTmp);
           % if max so far
          if tmp_e < e
               i = index;
               theta = thresh;
               p = pTmp;
               e = tmp_e;
               guess_y = tmp_y_guess;   
          end
        end
    end
end
   
%Try impossible values: check for decusion stump: return allways 1 / -1.
for pTmp=-1:2:1
  thresh = 3;
  [tmp_e, tmp_y_guess] = run_weak_classifier(X,y,w,index,thresh,pTmp);
   % if max so far
  if tmp_e < e
       i = index;
       theta = thresh;
       p = pTmp;
       e = tmp_e;
       guess_y = tmp_y_guess;   
  end
end

end


function [e, guess_y] = run_weak_classifier(X,y,w, i, theta,p)

% 1 for every sample that pass the decision stump. (equal to theta).
% -1 for the ones below it.
guess_y = 2*(X(i,:) == theta) - 1;

% consider p: if the decision stump is below threshold or above.
% also transpose it.
guess_y = guess_y' * p;

%check how many errors between the algo results and the real labels.
Errors = (guess_y ~= y);

% weight the errors with the w distrebution.
e = sum(Errors .* w);

end
