function [ algorithm ] = entangled_snp_algorithm( dim_radius )

algorithm.train = @train;
algorithm.classify = @classify;
algorithm.description = sprintf('Entangled-Snp, dim_radius = %d', dim_radius);
algorithm.params.dim_radius = dim_radius;

end

function [ model ] = train(params, train, extracted_train, snp_positions, missing, weights)

if nargin < 6
    weights = ones(size(extracted_train,3), 1) / size(extracted_train,3);
end
  
model.indices = zeros(1, length(missing));
model.perm = zeros(3, length(missing));

middle = (size(extracted_train,2) + 1 )/2;

if params.dim_radius == 0
    params.dim_radius = middle - 1;
end

extracted_train = extracted_train(:,middle-params.dim_radius:middle+params.dim_radius,:);

for i = 1:length(missing)
    [current_index, current_perm] = find_entangled_snp(squeeze(extracted_train(i,:,:)), weights);
    model.indices(i) = current_index;
    model.perm(:,i) = current_perm;
end

end

function [index, perm] = find_entangled_snp(extracted_train, weights)

all_perms = perms([0 1 2]);
indices = zeros(size(all_perms,1),1);
errors = zeros(size(all_perms,1),1);

middle = (size(extracted_train,1) + 1 )/2;
original_labels = extracted_train(middle,:);

extracted_train = changem(extracted_train, [3 4 5], [0 1 2]);
extracted_train(middle,:) = original_labels;

for i = 1:size(all_perms,1)
    [current_index, current_error] = find_entangled_snp_for_perm( ...
        changem(extracted_train, all_perms(i,:), [3 4 5]), ...
        weights);
    
    indices(i) = current_index;
    errors(i) = current_error;
end

[min_error, perm_index] = min(errors);
perm = all_perms(perm_index,:);
index =  indices(perm_index);

end

function [index, error] = find_entangled_snp_for_perm(extracted_train, weights) % TODO: use all data?

middle = (size(extracted_train,1) + 1 )/2;

ground_truth_matrix = repmat(extracted_train(middle,:),size(extracted_train,1),1);
weights_matrix = repmat(weights',size(extracted_train,1),1);

diff_matrix = ground_truth_matrix ~= extracted_train; % TODO: maybe difference? 1 is closer to 2 then 0 is
diff_matrix(middle,:) = ones(size(extracted_train,2),1);
diff_matrix = double(diff_matrix) .* weights_matrix;

error_rates = sum(diff_matrix, 2);

[error, index] = min(error_rates);

end

function [ ytest ] = classify(model, test, extracted_test, snp_positions, missing)

ytest = zeros(length(missing), size(test,2));

extracted_test = changem(extracted_test, [3 4 5], [0 1 2]);

for i = 1:length(missing)
    ytest(i,:) = changem(...
        squeeze(extracted_test(i,model.indices(i),:)),...
        model.perm(:,i),...
        [3 4 5]);
end

end