function [ algorithm ] = tree_algorithm(dim_radius, extra_param)

algorithm.train = @tree_train;
algorithm.classify = @tree_classify;
algorithm.description = sprintf('decission tree, options = %s', extra_param);
algorithm.params.dim_radius = dim_radius;
algorithm.params.extra_param = extra_param;

end

function [ model ] = tree_train(params, train, extracted_train, snp_positions, missing)

n = length(missing);

model.tree_models = cell(n,1);
model.dim_radius = params.dim_radius;

middle = (size(extracted_train,2) + 1 )/2;
extracted_train = extracted_train(:,middle-params.dim_radius:middle+params.dim_radius,:);

for i = 1:n
    model.tree_models{i} = train_single_tree_model(...
        squeeze(extracted_train(i,:,:)), ...
        params);
    if mod(i,50) == 0
        fprintf('Done: %d/%d \n', i, n);
    end
end

end

function [ model ] = train_single_tree_model(extracted_train_snp, params)

extra_param = params.extra_param;

middle = (size(extracted_train_snp,1) + 1 )/2;
train_labels = double(extracted_train_snp(middle, :))';
X = extracted_train_snp([1:middle-1,middle+1:end],:);

X = double(X);

model = ClassificationTree.fit(X',train_labels, 'AlgorithmForCategorical','Exact');

model = prune(model,'Level',2);

end

function [ ytest ] = tree_classify(model, test, extracted_test, snp_positions, missing)

snps = length(missing);
peoples = size(extracted_test,3);

ytest = zeros(snps, peoples);

middle = (size(extracted_test,2) + 1 )/2;
extracted_test = extracted_test(:,middle-model.dim_radius:middle+model.dim_radius,:);

for i = 1:snps   
    extracted_test_snp = squeeze(extracted_test(i,:,:));
    middle = (size(extracted_test_snp,1) + 1 )/2;
    X = extracted_test_snp([1:middle-1,middle+1:end],:);
    
    X = double(X);
    
    predictions = model.tree_models{i}.predict(X');
    ytest(i,:) = predictions';
end

end


