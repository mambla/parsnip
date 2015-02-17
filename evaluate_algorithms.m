function [ accuracy_vec, best_algorithm ] = evaluate_algorithms(...
    algorithms, ... % documentation...
    train, ...
    extracted_train, ...
    snp_positions, ...
    missing)

accuracy_vec = zeros(length(algorithms), 1);

for i = 1:length(algorithms)
    current_algorithm =  algorithms(i);
    
    fprintf('Testing %s algorithm \n', current_algorithm.description);
    
    accuracy_vec(i) = cross_validate(...
        current_algorithm, train, extracted_train, snp_positions, missing, 3);
    
    fprintf('Done testing %s algorithm. Accuracy: %f \n', ...
        current_algorithm.description, accuracy_vec(i));
    
    [~, current_best_algorithm_index] = max(accuracy_vec);
    current_best_algorithm = algorithms(current_best_algorithm_index);
    
    fprintf('Current best algorithm: %s, with accuracy: %f \n', ...
        current_best_algorithm.description, accuracy_vec(current_best_algorithm_index));
end

[~, best_algorithm_index] = max(accuracy_vec);
best_algorithm = algorithms(best_algorithm_index);

end

