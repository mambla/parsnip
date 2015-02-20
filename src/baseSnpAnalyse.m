% get only the missed snp's of training set.
% matrix of snp*samples
snps = squeeze(extracted_train(:,101,:));

numOfSNP=300;
x = zeros(numOfSNP,3);

for i = 1:numOfSNP
    for j=1:3
        x(i,j) = mean( snps(i,:) == j-1 );
    end
end

% x(i,j) == the statistic over train set, that the i special snp is going
% to be value <j-1>
%
% y is for each speial snp; what is the success rate when guessing the
% majoriy
y = max(x,[],1);