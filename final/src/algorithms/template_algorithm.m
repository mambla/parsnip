function [ algorithm ] = template_algorithm(const)
% a template for the struct used by all other algorithms
%
% every algorithm should contain a train and classify function 
% (according to the below signatures, a
% description, and parameters if needed
%
% this algorithm always guesses the constant given as an argument

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


