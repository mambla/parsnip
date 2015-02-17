function [ algorithm ] = template_algorithm( )

algorithm.train = @train;
algorithm.classify = @classify;
algorithm.description = sprintf('template');
algorithm.params = 0;

end

function [ model ] = train(params, train, extracted_train, snp_positions, missing)

model = 0;

end

function [ ytest] = classify(model, test, extracted_test, snp_positions, missing)

ytest = ones(length(missing), size(test,2));

end


