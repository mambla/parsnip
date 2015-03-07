function [ algorithm ] = template_algorithm(const)

algorithm.train = @template_train;
algorithm.classify = @template_classify;
algorithm.description = sprintf('template');
algorithm.params.const = const;

end

function [ model ] = template_train(params, train, extracted_train, snp_positions, missing)

model.const = params.const;

end

function [ ytest] = template_classify(model, test, extracted_test, snp_positions, missing)

ytest = ones(length(missing), size(test,2)) * model.const;

end


