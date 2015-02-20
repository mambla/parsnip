function [ Model ] = TrainAdaboost( X, y, target, nRounds)
%X's column is sample. c's row is feature.
% y labels vector (0/1/2)
% target=0/1/2. we would build classifier that detect this value or not
% this value
%nRounds - num of adaboost rounds

addpath('./fwdboosting/');
sPARAMS.Nrounds = nRounds;

switch target
    case 0
        y = changem(y,[1 -1 -1], [0 1 2]);
    case 1
        y = changem(y,[-1 1 -1], [0 1 2]);
    case 2
        y = changem(y,[-1 -1 1], [0 1 2]);
end

Model = CLSgentleBoost(X,y,sPARAMS);
end
