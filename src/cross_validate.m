function [ accuracy, accuracyV ] = cross_validate(...
    algorithm, ...  % documentation...
    train, ...
    extracted_train, ...
    snp_positions, ...
    missing, ...
    k)

% TODO: Split data and then try all algorithms. Data manipulation takes
% time because of size

% Take documentation from ex3 k_cross_validate...

% Initialize accuracy vector
accuracyV = zeros(length(missing),1);

train_size = size(train, 2);

% randomally permute data and calculate group splits
perm = randperm(train_size);
perm_train = train(:, perm);
perm_extracted_train = extracted_train(:, :, perm);

group_size = floor(train_size / k);
group_ranges = 1:group_size:train_size - k + 1;
group_ranges = [group_ranges train_size+1];

for i = 1:k
    % omit one group from the training set, and choose it as the test set
    train_set = perm_train;
    extracted_train_set = perm_extracted_train;
    train_set(:, group_ranges(i):group_ranges(i+1)-1) = [];
    extracted_train_set(:, :, group_ranges(i):group_ranges(i+1)-1) = [];
    
	test_set = perm_train(:, group_ranges(i):group_ranges(i+1)-1);
    extracted_test_set = perm_extracted_train(:, :, group_ranges(i):group_ranges(i+1)-1);
    
    % set SNPs in test set to -1  
    middle = (size(extracted_test_set,2) + 1 )/2;
    ground_truth = squeeze(extracted_test_set(:, middle, :));
    test_set(missing, :) = -1;
    extracted_test_set(:, middle, :) = -1;
    
    % train
    model = algorithm.train(...
        algorithm.params, train_set, extracted_train_set, snp_positions, missing);
    
    % classify
    ytest = algorithm.classify(...
        model, test_set, extracted_test_set, snp_positions, missing);
    
    % cumulate accuracy for average calculated at the end of the function
    tmp_accuracy =  mean(ytest ==  ground_truth(1:length(missing),:),2);
    accuracyV = accuracyV + tmp_accuracy;  
    fprintf('Cross Validation Done: %d/%d \t Accuracy:%f\n', i, k, mean(tmp_accuracy));

 
end

% calculate average accuracy
accuracyV = accuracyV / k;
accuracy = mean(accuracyV);

end

