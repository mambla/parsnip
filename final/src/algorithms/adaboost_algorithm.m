function [ algorithm ] = adaboost_algorithm(T, radius)

% <T> iterations adaboost upon the <radius> SNPs before and after
% the one we are looking for
% weak classifiers: does snp in index i has value x?
% T and radius are given as parameters

algorithm.train = @adaboost_train;
algorithm.classify = @adaboost_classify;
algorithm.description = sprintf('Adaboost, T = %d', T);
algorithm.params.T = T;
algorithm.params.radius = radius;

end

function [ model ] = adaboost_train(params, train, extracted_train, snp_positions, missing)

%initialize model
model.adaboost_0 = cell(length(missing),1);
model.adaboost_1 = cell(length(missing),1);
model.adaboost_2 = cell(length(missing),1);
model.T = params.T;
model.radius = params.radius;

for i = 1:length(missing)
    [ X, y ] = adaboost_init( extracted_train, i, params.radius );
    
    %0 detector
    model.adaboost_0{i} = TrainAdaboost2( X, y, 0, params.T);
    %1 detector
    model.adaboost_1{i} = TrainAdaboost2( X, y, 1, params.T);
    %2 detector
    model.adaboost_2{i} = TrainAdaboost2( X, y, 2, params.T);
    
    fprintf('Done: %d/%d \n', i, length(missing));
end

end

function [ ytest ] = adaboost_classify(model, test, extracted_test, snp_positions, missing)

ytest = zeros(length(missing), size(test,2));

%for each snp
for i = 1:length(missing)
    [ X, ~ ] = adaboost_init( extracted_test, i, model.radius );
    [ ~, r0 ] = adaboostC2( X, model.adaboost_0{i} );
    [ ~, r1 ] = adaboostC2( X, model.adaboost_1{i} );
    [ ~, r2 ] = adaboostC2( X, model.adaboost_2{i} );
    
    % in total: check who had maximal r. and give this label
    [~, max_r_indexes] = max([r0 r1 r2],[],2);
    ytest(i,:) = (max_r_indexes-1)';
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   sub routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ y_classifier_result, r ] = adaboostC2( X, Model )
% this run the classifier on the test samples.
% return as well r - how mutch the classifier is sure about it's answer

i_feature = Model.i;
thresh = Model.t;
p = Model.p;
alpha = Model.alpha;

% num of samples
n = size(X,2);

r = zeros(n,1);
y_classifier_result = zeros(n,1);

%for each sample
for k=1:n
    r(k) = sum(alpha .* ((2*(X(i_feature,k) == thresh)-1) .* p));
    y_classifier_result(k) = sign(r(k));
end

end




function [ Model ] = TrainAdaboost2( X, y, targetClass, nRounds)
%X's column is sample. c's row is feature.
% y labels vector (0/1/2)
% target=0/1/2. we would build classifier that detect this value or not
% this value
%nRounds - num of adaboost rounds

sPARAMS.T = nRounds;

switch targetClass
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


function [ X, y ] = adaboost_init( extracted_train, snp_ind, r )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    r = r - 1;
    tmp = squeeze(extracted_train(snp_ind,:,:));
    X = tmp([100-r:100,102:102+r],:);
    y = tmp(101,:);
    y = y';
    %clearvars tmp;
    
    X = double(X);
    y = double(y);
end



