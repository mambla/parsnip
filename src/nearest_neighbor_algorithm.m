function [ algorithm ] = nearest_neighbor_algorithm( k, sigma)

%TODO: maybe weighted kNN, by position or index (e^-abs(Dist)), 0 weight
%for different chromosomes
%TODO: try different distance metrics
%TODO: use matlab knn
%TODO: different types of weight vectors (rect for instance)
%TODO: different sigma per snp

algorithm.train = @train;
algorithm.classify = @classify;
algorithm.description = sprintf('%d-NN, sigma=%d', k, sigma);
algorithm.params.k = k;
algorithm.params.sigma = sigma;

end

function [ model ] = train(params, train, extracted_train, snp_positions, missing)

model.extracted_train = extracted_train;
model.k = params.k;
model.sigma = params.sigma;

end

function [ ytest] = classify(model, test, extracted_test, snp_positions, missing)

ytest = zeros(length(missing), size(test,2));
weights = fspecial('gaussian', [201 1], model.sigma);

for i = 1:length(missing)
    ytest(i,:) = classify_single_snp(...
        squeeze(model.extracted_train(i,:,:)), ...
        squeeze(extracted_test(i,:,:)), ...
        model.k, ...
        weights);
end

end

function [ yrow ] = classify_single_snp(extracted_train_snp, extracted_test_snp, k, weights)

train_labels = extracted_train_snp(101, :);
extracted_train_snp(101, :) = -1;

extracted_train_snp = double(extracted_train_snp).*repmat(weights, 1, size(extracted_train_snp,2));
extracted_test_snp = double(extracted_test_snp).*repmat(weights, 1, size(extracted_test_snp,2));

dist_matrix = pdist2(double(extracted_train_snp)', double(extracted_test_snp)');

% sort every column in the dist matrix and extract the sorting permutation
[~, sorted_index_matrix] = sort(dist_matrix, 1, 'ascend');

% extract the labels for the k closest neighbor for every row of the test
% data by choosing the labels from class_T according to the first k values
% of the sorting permutation
sorted_label_matrix = train_labels(sorted_index_matrix(1:k,:)');

% calculate the plurality for each row of k closes labels. If there is a
% tie, the row contains a cell array with all possible values
%[~,~,yrow_with_ties] = mode(sorted_label_matrix,2);
if k > 1
    yrow = mode(sorted_label_matrix,2);
else
    yrow = sorted_label_matrix;
end

% random tiebreak: choose a random value from every cell array in the
% plurality result from before
% for i = 1:size(yrow_with_ties,1)
%    current_tiebreak = yrow_with_ties{i};
%    yrow(i,1) = current_tiebreak(randi(size(current_tiebreak,1)));
% end


end



