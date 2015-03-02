top_features = 10;
radius = 40;
%svm_options = '-t 2 -g 0.01 -c 100 -q';
svm_options = '-t 2 -g 0.01 -c 100 -q';

svm2 = svm_algorithm2(svm_options,radius,top_features);

K=3;
%save time, run all the process only on snp 1 to <snp_to_test>
snp_to_test = 300;
%the adaboost doesn't use the all train or all test. only the extracted
%data. I passed dummy matrixes to save time.
[~, svmRate] = cross_validate(svm2,train(1,:),extracted_train,[],missing(1:snp_to_test),K);

%Error is mean error for all snp
%errorV is the mean error for each separate snp tested