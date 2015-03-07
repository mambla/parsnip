function [ ytest, chosen_algorithms, accuracy_vector, accuracy ] = go()
% chooses the best algorithm for each snp (using 6-fold cross validation)
% and estimates the missing SNPs in the test set using the chosen
% algorithms
%
% this function returns and saves ytest (in ytest.mat) , the chosen_algorithms, an
% accuracy_vector that contains the cross validation accuracy of the
% algorithm chosen for each SNP, and the average accuracy (all in
% chosen_algorithms.mat)

% loads data and adds paths
setup();

% if algorithms were previously chosen per SNP, load the previously chosen
% algorithms, otherwise, find the best algorithm per SNP
if exist('../results/chosen_algorithms.mat','file')
    load chosen_algorithms.mat
else
    [chosen_algorithms, accuracy_vector, accuracy] = choose_algorithm_per_snp(...
        train, extracted_train, missing);
    
    % save for next time
    save('../results/chosen_algorithms.mat','chosen_algorithms','accuracy_vector','accuracy');
end

% initialize empty result
ytest = zeros(length(missing), size(test,2));

for i = 1:length(missing)
    % train and then classify using previously chosen algorithms
    current_algorithm = chosen_algorithms(i);
    
    current_model = current_algorithm.train(...
        current_algorithm.params, train, extracted_train(i,:,:), 0, missing(i));
    
    ytest(i,:) = current_algorithm.classify(...
        current_model, test, extracted_test(i,:,:), 0, missing(i));
end

% save result
save('../results/ytest.mat','ytest');

end

function [chosen_algorithms, accuracy_vector, accuracy] = choose_algorithm_per_snp(...
    train, extracted_train, missing)

% get algorithm options
algorithm_options = generate_algorithm_options();

% create matrix that will eventually contain cross validation accuracy for
% each SNP with each algorithm
accuracy_matrix = zeros(length(missing), length(algorithm_options));

% calculate cross validation accuracy for each algorithm on all SNPs
for i = 1:length(algorithm_options)
    fprintf('Testing algorithm: %s \n', current_algorithm.description);
    
    % snp positions are not used in any algorithm (the 0 argument)
    [~, current_accuracy_vector] = cross_validate(...
        algorithm_options(i), train, extracted_train, 0, missing, 6);
    
    accuracy_matrix(:,i) = current_accuracy_vector;
    
    fprintf('Done testing algorithm: %s. Accuracy: %f \n', ...
        current_algorithm.description, accuracy_vec(i));
end

% find max values / indexes of algorithms that acheived the best accuracy
% for each SNP and extract chosen algorithms
[accuracy_vector, algo_index_per_snp] = max(accuracy_matrix, [], 2);
chosen_algorithms = algorithm_options(algo_index_per_snp);

% calculate mean accuracy
accuracy = mean(accuracy_vector);

end

function [algorithm_options] = generate_algorithm_options()

% initialize a vector of the algorithms (with different paramters) from which we wish to choose
% the algorithm parameters were chosen as parameters that previously
% produced good results

algorithm_options = [...
    svm_algorithm('-t 2 -g 0.007 -c 100 -q',9), ...
    entangled_snp_algorithm(20), ...
    nearest_neighbor_algorithm(3,3, ones(201,1))];
    
end