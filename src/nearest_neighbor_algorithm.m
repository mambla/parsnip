function [ algorithm ] = nearest_neighbor_algorithm( k, sigma)

%TODO: maybe weighted kNN, by position or index (e^-abs(Dist)), 0 weight
%for different chromosomes
%TODO: try different distance metrics
%TODO: use matlab knn
%TODO: different types of weight vectors (rect for instance)
%TODO: different sigma per snp
%TODO: different algorithm per snp?

algorithm.train = @train;
algorithm.classify = @classify;
algorithm.description = sprintf('%d-NN, sigma=%d', k, sigma);
algorithm.params.k = k;
algorithm.params.sigma = sigma;

end

function [ model ] = train(params, train, extracted_train, snp_positions, missing)

model.k = params.k;
model.weights = fspecial('gaussian', [size(extracted_train,2) 1], params.sigma);
model.knn_mdls = cell(length(missing),1);

for i = 1:length(missing)
   model.knn_mdls{i} = train_single_snp_model(...
        squeeze(extracted_train(i,:,:)), ...
        model.k, ...
        model.weights);
end

end

function [ ytest ] = classify(model, test, extracted_test, snp_positions, missing)

ytest = zeros(length(missing), size(test,2));

for i = 1:length(missing)
    extracted_test_snp = squeeze(extracted_test(i,:,:));
    extracted_test_snp = double(extracted_test_snp).*repmat(model.weights, 1, size(extracted_test_snp,2));

    ytest(i,:) = predict(model.knn_mdls{i}, extracted_test_snp')';
end

end

function [ mdl ] = train_single_snp_model(extracted_train_snp, k, weights)

middle = (size(extracted_train_snp,1) + 1 )/2;
train_labels = extracted_train_snp(middle, :);
extracted_train_snp(middle, :) = -1;

extracted_train_snp = double(extracted_train_snp).*repmat(weights, 1, size(extracted_train_snp,2));

mdl = fitcknn(extracted_train_snp', train_labels','NumNeighbors',k);%, 'Distance', 'correlation');

end



