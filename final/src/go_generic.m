function [ ytest ] = go_generic(algorithm)
% estimates ytest (saves and returns) using the given algorithm

% loads data and adds paths
setup();

fprintf('Using algorithm: %s\n', algorithm.description);
    
% train algorithm, snp positions aren't used in any algorithm (the 0 argument)
model = algorithm.train(algorithm.params, train, extracted_train, 0, missing);   

% estimate missing SNPs
ytest = algorithm.classify(model, test, extracted_test, 0, missing);
    
end