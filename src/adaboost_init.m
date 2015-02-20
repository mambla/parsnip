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


