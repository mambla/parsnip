function [ algorithm ] = svm_algorithm2(svm_options, dim_radius, top_features)

algorithm.train = @train;
algorithm.classify = @classify;
algorithm.description = sprintf('svm, options = %s', svm_options);
algorithm.params.svm_options = svm_options;
algorithm.params.dim_radius = dim_radius;
algorithm.params.top_features = top_features;

end

% successful params: svm_algorithm('-t 2 -g 0.007 -c 100 -q',20)
function [ model ] = train(params, train, extracted_train, snp_positions, missing)

n = length(missing);

model.svm_models = cell(n,1);
model.features = cell(n,1);
model.dim_radius = params.dim_radius;

middle = (size(extracted_train,2) + 1 )/2;
extracted_train = extracted_train(:,middle-params.dim_radius:middle+params.dim_radius,:);

for i = 1:n
    [model.svm_models{i}, model.features{i}]  = train_single_svm_model(...
        squeeze(extracted_train(i,:,:)), ...
        params);
    if mod(i,50) == 0
        fprintf('Done: %d/%d \n', i, n);
    end
end

end

function [ model, features] = train_single_svm_model(extracted_train_snp, params)

top_features = params.top_features;
svm_options = params.svm_options;

middle = (size(extracted_train_snp,1) + 1 )/2;
train_labels = double(extracted_train_snp(middle, :))';
X = extracted_train_snp([1:middle-1,middle+1:end],:);
X = double(X);
% find meaningful features

% check all permutations of [0,1,2]
XX = [X; changem(X,[0 2 1], [0 1 2]); changem(X,[1 0 2], [0 1 2]); changem(X,[1 2 0], [0 1 2]); changem(X,[2 0 1], [0 1 2]); changem(X,[2 1 0], [0 1 2])];

% calculate correlation of labels with the permutated near (dist< r) snps
cor = zeros(size(XX,1), 1);
for i=1:size(cor,1)
    cor(i) = mean(train_labels' == XX(i,:));
end

%keep only meaningful features
[~,ind] = sort(cor,'descend');
%[~,ind] = sort(cor,'ascend'); %???

features = ind(1:top_features);
Train_features = XX(features,:);

model = svmtrain(train_labels, Train_features', svm_options);

end

function [ ytest ] = classify(model, test, extracted_test, snp_positions, missing)

snps = length(missing);
peoples = size(extracted_test,3);

ytest = zeros(snps, peoples);

middle = (size(extracted_test,2) + 1 )/2;
extracted_test = extracted_test(:,middle-model.dim_radius:middle+model.dim_radius,:);

for i = 1:snps   
    extracted_test_snp = squeeze(extracted_test(i,:,:));
    middle = (size(extracted_test_snp,1) + 1 )/2;
    X = extracted_test_snp([1:middle-1,middle+1:end],:);
    XX = [X; changem(X,[0 2 1], [0 1 2]); changem(X,[1 0 2], [0 1 2]); changem(X,[1 2 0], [0 1 2]); changem(X,[2 0 1], [0 1 2]); changem(X,[2 1 0], [0 1 2])];
    Test_features = XX(model.features{i},:);
    
    ytest(i,:) = svmpredict(ytest(i,:)', double(Test_features'), model.svm_models{i}, '-q')';
end

end


