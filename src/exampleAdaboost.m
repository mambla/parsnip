T = 9;
radius = 30;
ada = adaboost_algorithm(T,radius);

K=6;
%save time, run all the process only on snp 1 to <snp_to_test>
snp_to_test = 300;
%the adaboost doesn't use the all train or all test. only the extracted
%data. I passed dummy matrixes to save time.
[error, errorV] = cross_validate(ada,train(1,:),extracted_train,[],missing(1:snp_to_test),K);

%Error is mean error for all snp
%errorV is the mean error for each separate snp tested