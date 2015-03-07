function [ algorithm ] = entangled_snp_algorithm( dim_radius )
% This algorithm searches for the most "entangled" SNP to each missing SNP
% The most entangled SNP is found to be the SNP (or any of the 27 substitutions 
% of the SNP, [0 1 2] -> [? ? ?],  within the dim_radius before and after the 
% missing SNP whose percentage of different values is lowest. Once found,
% classifying is done by simply choosing the value of the entangled SNP
% after substitution

algorithm.train = @entangled_snp_train;
algorithm.classify = @entangled_snp_classify;
algorithm.description = sprintf('entangled_snp, dim_radius = %d', dim_radius);
algorithm.params.dim_radius = dim_radius;

end

function [ model ] = entangled_snp_train(params, train, extracted_train, snp_positions, missing)

model.indices = zeros(1, length(missing));
model.perm = zeros(3, length(missing));

% find middle of data, the missing SNP
middle = (size(extracted_train,2) + 1 )/2;

% if dim radius is 0, use all data in extracted_train
if params.dim_radius == 0
    params.dim_radius = middle - 1;
end

model.dim_radius = params.dim_radius;

% extract dim_radius SNPs before and after the middle SNP
extracted_train = extracted_train(:,middle-params.dim_radius:middle+params.dim_radius,:);

% find entangled SNP for each SNP
for i = 1:length(missing)
    [current_index, current_perm] = find_entangled_snp(squeeze(extracted_train(i,:,:)));
    model.indices(i) = current_index;
    model.perm(:,i) = current_perm;
end

end

function [index, perm] = find_entangled_snp(extracted_train)

% create all transformations: [0 1 2] -> [? ? ?]
all_perms = unique(combnk([0:2 0:2 0:2],3),'rows');
indices = zeros(size(all_perms,1),1);
errors = zeros(size(all_perms,1),1);

middle = (size(extracted_train,1) + 1 )/2;
original_labels = extracted_train(middle,:);

% transform [0 1 2] to [3 4 5] except for the SNP we are trying to predict
extracted_train = changem(extracted_train, [3 4 5], [0 1 2]);
extracted_train(middle,:) = original_labels;

% find most entangled SNP for each substitution
for i = 1:size(all_perms,1)
    % find entangled SNP after substitution
    [current_index, current_error] = find_entangled_snp_for_perm( ...
        changem(extracted_train, all_perms(i,:), [3 4 5]));
    
    indices(i) = current_index;
    errors(i) = current_error;
end

% choose best substitution (lowest error on training set)
[~, perm_index] = min(errors);
perm = all_perms(perm_index,:);
index =  indices(perm_index);

end

function [index, error] = find_entangled_snp_for_perm(extracted_train)

middle = (size(extracted_train,1) + 1 )/2;

% replicate labels
ground_truth_matrix = repmat(extracted_train(middle,:),size(extracted_train,1),1);

% compare labels and train values
diff_matrix = ground_truth_matrix ~= extracted_train; 

% avoid choosing the labels of the missing SNP as the most entangled (with itself)
% by placing one at values of diff matrix where the missing SNP is
diff_matrix(middle,:) = ones(size(extracted_train,2),1);

% sum of the matrix is the number of errors for each SNP
error_rates = sum(diff_matrix, 2);

% choose the one with the least amount of errors
[error, index] = min(error_rates);

end

function [ ytest ] = entangled_snp_classify(model, test, extracted_test, snp_positions, missing)

ytest = zeros(length(missing), size(test,2));

% choose the most entangled SNP after substitution for each missing SNP 

middle = (size(extracted_test,2) + 1 )/2;
extracted_test = extracted_test(:,middle-model.dim_radius:middle+model.dim_radius,:);
extracted_test = changem(extracted_test, [3 4 5], [0 1 2]);

for i = 1:length(missing)
    ytest(i,:) = changem(...
        squeeze(extracted_test(i,model.indices(i),:)),...
        model.perm(:,i),...
        [3 4 5]);
end

end