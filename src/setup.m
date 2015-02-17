addpath('.');
addpath('../data');
addpath('../bin');
addpath('../bin/windows');

load('dataforproject.mat');

if ~exist('../data/train.mat','file')
    train = dlmread('../data/train.txt');
    test = dlmread('../data/test.txt');
    save('../data/train.mat','train');
    save('../data/test.mat','test');
else
    load('../data/train.mat','train');
    load('../data/test.mat','test');
end
