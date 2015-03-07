function [ algorithm ] = nearest_neighbor_algorithm( k, sigma, distance_metric)
% This algorithm uses KNN upon extracted train values. We found that better
% results were acheived by multiplying the values of the SNPs by a guassian
% weights vector (with the given sigma as the variance) rather then using 
% only a limited amount of values surrounding the missing SNP. The k and
% distance metric of the KNN are also given as arguments

algorithm.train = @nearest_neighbor_train;
algorithm.classify = @nearest_neighbor_classify;
algorithm.description = sprintf('%d-NN, sigma=%d, distance_metric=%s', k, sigma, distance_metric);
algorithm.params.k = k;
algorithm.params.sigma = sigma;
algorithm.params.distance_metric = distance_metric;

end

function [ model ] = nearest_neighbor_train(params, train, extracted_train, snp_positions, missing)

model.k = params.k;

%calculate weights vector
model.weights = fspecial('gaussian', [size(extracted_train,2) 1], params.sigma);

% initialize model cell array
model.knn_mdls = cell(length(missing),1);

% train KNN model for each SNP
for i = 1:length(missing)
   model.knn_mdls{i} = train_single_snp_model(...
        squeeze(extracted_train(i,:,:)), ...
        model.k, ...
        model.weights,...
        params.distance_metric);
    
end

end

function [ ytest ] = nearest_neighbor_classify(model, test, extracted_test, snp_positions, missing)

ytest = zeros(length(missing), size(test,2));

% multiply the test data of each SNP by the weights vector and then classify using the
% knn model
for i = 1:length(missing)
    extracted_test_snp = squeeze(extracted_test(i,:,:));
    extracted_test_snp = double(extracted_test_snp).*repmat(model.weights, 1, size(extracted_test_snp,2));

    ytest(i,:) = predict(model.knn_mdls{i}, extracted_test_snp')';
end

end

function [ mdl ] = train_single_snp_model(extracted_train_snp, k, weights, distance_metric)

% extract labels
middle = (size(extracted_train_snp,1) + 1 )/2;
train_labels = extracted_train_snp(middle, :);
extracted_train_snp(middle, :) = -1;

% multiply by weights
extracted_train_snp = double(extracted_train_snp).*repmat(weights, 1, size(extracted_train_snp,2));

% train knn model using matlab implementation
mdl = fitcknn(extracted_train_snp', train_labels','NumNeighbors',k, 'Distance', distance_metric);

end



