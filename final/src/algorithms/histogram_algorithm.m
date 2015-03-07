function [ algorithm ] = histogram_algorithm(svm_options, width, slide_interval)
% This algorithm uses SVM trained using a feature vector that contains
% the concatenation of histograms of different sequences of SNPs withing the extracted data. 
% The sequences the histrogram is calculated upon is defined by the arguments width
% (amount SNPs to calculate the histogram upon) and slide interval (the
% interval of the sliding window of width, width, that the histrograms are
% calculated upon)

algorithm.train = @histogram_train;
algorithm.classify = @histogram_classify;
algorithm.description = sprintf('histogram, width = %d, slide_interval = %d, svm_options = %s', width, slide_interval, svm_options);
algorithm.params.svm_options = svm_options;
algorithm.params.width = width;
algorithm.params.slide_interval = slide_interval;

end

function [ model ] = histogram_train(params, train, extracted_train, snp_positions, missing)

% train svm model for each SNP
model.svm_models = cell(length(missing),1);
model.width = params.width;
model.slide_interval = params.slide_interval;

for i = 1:length(missing)
    model.svm_models{i} = train_single_svm_model(...
        model, ...
        squeeze(extracted_train(i,:,:)), ...
        params.svm_options);
    
    if mod(i,50) == 0
        fprintf('Done: %d/%d \n', i, length(missing));
    end
end

end

function [ hd ] = calculate_histogram_data(model, extracted_snp_data)
% calculate histogram sequence indices
indices = -model.width/2:model.slide_interval:size(extracted_snp_data,1)/2-model.width;
hd = zeros(length(model.width) * length(indices) * 2 * 3 ,size(extracted_snp_data,2));
middle = size(extracted_snp_data,1)/2;
current = 1;

% Calculate histograms in a symetric manner. First a histrogram is
% calculated upon the values in the middle of the extraced train data, and
% then we slide, slide_interval SNPs in each direction and calculate the
% next histogram
for i = indices
    % slide up
    current_data = extracted_snp_data(middle+i+1:middle+i+model.width,:);
    hd(current:current + 2,:) = hist(single(current_data),3);
    current = current + 3;
    
    % slide down
    current_data = extracted_snp_data(middle-i-model.width+1:middle-i,:);
    hd(current:current + 2,:) = hist(single(current_data),3);
    current = current + 3;
end

end

function [ model ] = train_single_svm_model(model, extracted_train_snp, svm_options)
% exactly like the SVM algorithm, but the histogram data is used
middle = (size(extracted_train_snp,1) + 1 )/2;
train_labels = extracted_train_snp(middle, :);
extracted_train_snp(middle, :) = [];
hd = calculate_histogram_data(model, extracted_train_snp);

model = svmtrain(double(train_labels)', double(hd)', svm_options);

end

function [ ytest ] = histogram_classify(model, test, extracted_test, snp_positions, missing)
% exactly like the SVM algorithm, but the histogram data is used
ytest = zeros(length(missing), size(test,2));

middle = (size(extracted_test,2) + 1 )/2;
extracted_test(:,middle,:) = [];

for i = 1:length(missing)
    extracted_test_snp = squeeze(extracted_test(i,:,:));
    hd = calculate_histogram_data(model, extracted_test_snp);
    ytest(i,:) = svmpredict(ytest(i,:)', double(hd'), model.svm_models{i}, '-q')';
end

end


