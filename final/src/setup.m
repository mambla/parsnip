% Add folders to search path
addpath 'algorithms';
addpath 'diagnostics';

addpath('.');
addpath '../data';
addpath '../results';
addpath '../bin';
addpath '../bin/windows';

if ispc()
    addpath '../bin/libsvm/windows';
else
    addpath '../bin/libsvm/linux';
end

% load data
load dataforproject.mat

if ~exist('../data/train.mat','file')
    train = dlmread('../data/train.txt');
    test = dlmread('../data/test.txt');
    save('../data/train.mat','train');
    save('../data/test.mat','test');
else
    load('../data/train.mat','train');
    load('../data/test.mat','test');
end