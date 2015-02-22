function [ algorithm ] = alt_adaboost_algorithm( t, dim_radius )

algorithm.train = @train;
algorithm.classify = @classify;
algorithm.description = sprintf('Alt-Adaboost, t = %d', t);
algorithm.params.t = t;
algorithm.params.dim_radius = dim_radius;

end

function [ model ] = train(params, train, extracted_train, snp_positions, missing)

model.t = params.t;
model.dim_radius = params.dim_radius;

end

function [ ytest ] = classify(model, test, extracted_test, snp_positions, missing)

end