function [ acc ] = preSVM(extracted_train, snp_ind)

%snp_ind = 1;
data = extracted_train;
r = 15;
top_features = 10;
svm_op = '-t 2 -g 0.007 -c 100 -q';

r = r - 1;
tmp = squeeze(extracted_train(snp_ind,:,:));
X = tmp([100-r:100,102:102+r],:);
y = tmp(101,:);
y = y';
%clearvars tmp;

X = double(X);
y = double(y);

% find meaningful features

% check all permutations of [0,1,2]
XX = [X; changem(X,[0 2 1], [0 1 2]); changem(X,[1 0 2], [0 1 2]); changem(X,[1 2 0], [0 1 2]); changem(X,[2 0 1], [0 1 2]); changem(X,[2 1 0], [0 1 2])];

% calculate correlation of labels with the permutated near (dist< r) snps
cor = zeros(size(XX,1), 1);
for i=1:size(cor,1)
    cor(i) = mean(y' == XX(i,:));
end

%keep only meaningful features
[~,ind] = sort(cor,'descend');
XXX = XX(ind(1:top_features),:);

%SVM

Train_features = XXX(:,1:400);
Test_features = XXX(:,401:600);
Train_y = y(1:400);
Test_y = y(401:600);

model = svmtrain(Train_y, Train_features', svm_op);
[yy, acc, ~] = svmpredict(Test_y, Test_features', model, '-q');

acc= acc(1);
end