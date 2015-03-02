function [ algorithm ] = svm_algorithm3(svm_options, dim_radius)

algorithm.train = @train;
algorithm.classify = @classify;
algorithm.description = sprintf('svm, options = %s, dim_radius = %d', svm_options, dim_radius);
algorithm.params.svm_options = svm_options;
algorithm.params.dim_radius = dim_radius;

end

% successful params: svm_algorithm('-t 2 -g 0.007 -c 100 -q',20)
function [ model ] = train(params, train, extracted_train, snp_positions, missing)

model.svm_models = cell(length(missing),1);
model.dim_radius = params.dim_radius;

middle = (size(extracted_train,2) + 1 )/2;

extracted_train = extracted_train(:,middle-params.dim_radius:middle+params.dim_radius,:);

for i = 1:length(missing)
    model.svm_models{i} = train_single_svm_model(...
        squeeze(extracted_train(i,:,:)), ...
        params.svm_options);
    
    if mod(i,10) == 0
        fprintf('Done: %d/%d \n', i, length(missing));
    end
end

end

function [ model ] = train_single_svm_model(extracted_train_snp, svm_options)

middle = (size(extracted_train_snp,1) + 1 )/2;
train_labels = extracted_train_snp(middle, :);
extracted_train_snp(middle, :) = -1;

value_counts = [sum(extracted_train_snp==0);sum(extracted_train_snp==1);sum(extracted_train_snp==2)];

model = svmtrain(double(train_labels)', double(value_counts)', svm_options);

end

function [ ytest ] = classify(model, test, extracted_test, snp_positions, missing)

ytest = zeros(length(missing), size(test,2));


middle = (size(extracted_test,2) + 1 )/2;
extracted_test = extracted_test(:,middle-model.dim_radius:middle+model.dim_radius,:);

for i = 1:length(missing)
    extracted_test_snp = squeeze(extracted_test(i,:,:));
    
    value_counts = [sum(extracted_test_snp==0);sum(extracted_test_snp==1);sum(extracted_test_snp==2)];

    ytest(i,:) = svmpredict(ytest(i,:)', double(value_counts'), model.svm_models{i}, '-q')';
end

end


