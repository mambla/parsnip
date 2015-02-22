function [ algorithm ] = adaboost_algorithm(T, radius)

algorithm.train = @train;
algorithm.classify = @classify;
algorithm.description = sprintf('Adaboost, T = %d', T);
algorithm.params.T = T;
algorithm.params.radius = radius;

end

function [ model ] = train(params, train, extracted_train, snp_positions, missing)

model.adaboost_0 = cell(length(missing),1);
model.adaboost_1 = cell(length(missing),1);
model.adaboost_2 = cell(length(missing),1);
model.T = params.T;
model.radius = params.radius;

for i = 1:length(missing)
    [ X, y ] = adaboost_init( extracted_train, i, params.radius );
    
    %0 detector
    model.adaboost_0{i} = TrainAdaboost2( X, y, 0, params.T);
    %1 detector
    model.adaboost_1{i} = TrainAdaboost2( X, y, 1, params.T);
    %2 detector
    model.adaboost_2{i} = TrainAdaboost2( X, y, 2, params.T);
    
    fprintf('Done: %d/%d \n', i, length(missing));
end

end

function [ ytest ] = classify(model, test, extracted_test, snp_positions, missing)

ytest = zeros(length(missing), size(test,2));

%for each snp
for i = 1:length(missing)
    [ X, ~ ] = adaboost_init( extracted_test, i, model.radius );
    [ ~, r0 ] = adaboostC2( X, model.adaboost_0{i} );
    [ ~, r1 ] = adaboostC2( X, model.adaboost_1{i} );
    [ ~, r2 ] = adaboostC2( X, model.adaboost_2{i} );
    
    % in total: check who had maximal r. and give this label
    [~, max_r_indexes] = max([r0 r1 r2],[],2);
    ytest(i,:) = (max_r_indexes-1)';
end

end


