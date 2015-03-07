function [ algorithm ] = svm_algorithm(svm_options, dim_radius)

% algorithm that uses multiclass svm in order to estimate missing SNPs
% the svm_options and dim_radius (SNPs before and after to use in the
% feature vector we train the svm on) are given as arguments

algorithm.train = @svm_train;
algorithm.classify = @svm_classify;
algorithm.description = sprintf('svm, options = %s', svm_options);
algorithm.params.svm_options = svm_options;
algorithm.params.dim_radius = dim_radius;

end

function [ model ] = svm_train(params, train, extracted_train, snp_positions, missing)

% initialize model
model.svm_models = cell(length(missing),1);
model.dim_radius = params.dim_radius;

% calculate middle of data (the missing SNP)
middle = (size(extracted_train,2) + 1 )/2;

% extract dim_radius SNPs before and after the middle SNP
extracted_train = extracted_train(:,middle-params.dim_radius:middle+params.dim_radius,:);

fprintf('Starting to train %d svm classifiers \n', length(missing));

% train SVM classifier for each SNP separately
for i = 1:length(missing)
    model.svm_models{i} = train_single_svm_model(...
        squeeze(extracted_train(i,:,:)), ...
        params.svm_options);
    
    if mod(i,50) == 0
        fprintf('Done: %d/%d \n', i, length(missing));
    end
end

end

function [ model ] = train_single_svm_model(extracted_train_snp, svm_options)

% extract train labels and then remove them from data
middle = (size(extracted_train_snp,1) + 1 )/2;
train_labels = extracted_train_snp(middle, :);
extracted_train_snp(middle, :) = [];

% train svm
model = svmtrain(double(train_labels)', double(extracted_train_snp)', svm_options);

end

function [ ytest ] = svm_classify(model, test, extracted_test, snp_positions, missing)

ytest = zeros(length(missing), size(test,2));

% extract dim_radius SNPs before and after the middle SNP
middle = (size(extracted_test,2) + 1 )/2;
extracted_test = extracted_test(:,middle-model.dim_radius:middle+model.dim_radius,:);

% remove missing SNP
middle = (size(extracted_test,2) + 1 )/2;
extracted_test(:,middle,:) = [];

% classify each SNP using respective model calculated in the train phase
for i = 1:length(missing)
    extracted_test_snp = squeeze(extracted_test(i,:,:));
    ytest(i,:) = svmpredict(ytest(i,:)', double(extracted_test_snp'), model.svm_models{i}, '-q')';
end

end


