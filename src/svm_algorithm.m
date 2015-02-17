function [ algorithm ] = svm_algorithm(svm_options)

algorithm.train = @train;
algorithm.classify = @classify;
algorithm.description = sprintf('svm, options = %s', svm_options);
algorithm.params.svm_options = svm_options;

end

% successful params: [ accuracy_vec, best_algorithm ] = evaluate_algorithms(svm_algorithm('-t 2 -g 0.007 -c 10 -q'), train, extracted_train(:,91:111,:), 0, missing);
function [ model ] = train(params, train, extracted_train, snp_positions, missing)

model.svm_models = cell(length(missing),1);

for i = 1:length(missing)
    model.svm_models{i} = train_single_svm_model(...
        squeeze(extracted_train(i,:,:)), ...
        params.svm_options);
    
    fprintf('Done: %d/%d \n', i, length(missing));
end

end

function [ model ] = train_single_svm_model(extracted_train_snp, svm_options)

middle = (size(extracted_train_snp,1) + 1 )/2;
train_labels = extracted_train_snp(middle, :);
extracted_train_snp(middle, :) = -1;

model = svmtrain(double(train_labels)', double(extracted_train_snp)', svm_options);

end

function [ ytest ] = classify(model, test, extracted_test, snp_positions, missing)

ytest = zeros(length(missing), size(test,2));

for i = 1:length(missing)
    extracted_test_snp = squeeze(extracted_test(i,:,:));
    ytest(i,:) = svmpredict(ytest(i,:)', double(extracted_test_snp'), model.svm_models{i}, '-q')';
end

end


