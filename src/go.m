function [ ytest, chosen_algorithms, accuracy_vector, accuracy ] = go()

load dataforproject.mat
load train.mat;
load test.mat;
    
if exist('chosen_algorithms.mat','file')
    load chosen_algorithms.mat
else
    [chosen_algorithms, accuracy_vector, accuracy] = choose_algorithm_per_snp(...
        train, extracted_train, missing);
    
    save chosen_algorithms.mat chosen_algorithms accuracy_vector accuracy
end

ytest = zeros(length(missing), size(test,2));

for i = 1:length(missing)
    current_algorithm = chosen_algorithms(i);
    
    current_model = current_algorithm.train(...
        current_algorithm.params, train, extracted_train(i,:,:), 0, missing(i));
    
    ytest(i,:) = current_algorithm.classify(...
        current_model, test, extracted_test(i,:,:), 0, missing(i));
end

save ytest.mat ytest

end

function [chosen_algorithms, accuracy_vector, accuracy] = choose_algorithm_per_snp(...
    train, extracted_train, missing)

algorithm_options = generate_algorithm_options();
accuracy_matrix = zeros(length(missing), length(algorithm_options));

for i = 1:length(algorithm_options)
    [~, current_accuracy_vector] = cross_validate(...
        algorithm_options(i), train, extracted_train, 0, missing, 6);
    
    accuracy_matrix(:,i) = current_accuracy_vector;
end

[accuracy_vector, algo_index_per_snp] = max(accuracy_matrix, [], 2);
accuracy = mean(accuracy_vector);
chosen_algorithms = algorithm_options(algo_index_per_snp);

end

function [algorithm_options] = generate_algorithm_options()

algorithm_options = [...
    svm_algorithm('-t 2 -g 0.007 -c 100 -q',9), ...
    entangled_snp_algorithm(20), ...
    nearest_neighbor_algorithm(3,3, ones(201,1))];
    
end