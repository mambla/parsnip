snps = squeeze(extracted_train(:,101,:));
numOfSNP=300;
Total = 165229;

%find best correlation between the i snmp and the j
% pdist -> M(i,j) := the distance between #i snp and #j snp in training set

%D = squareform(pdist(snp));
%D2 = D + diag(1000*ones(1,300));

%D = zeros(numOfSNP, numOfSNP);
%for i=1:numOfSNP
%    for j=1:numOfSNP
%        D(i,j) = mean(snps(i,:)==snps(j,:));
%    end
%    D(i,i)=0; % I dont care that vi is ident to vi
%end

D = zeros(numOfSNP, Total);
for i=1:numOfSNP
   for j=1:Total
       D(i,j) = mean(snps(i,:)==train(j,:));  
   end
   % dont want someone to be correlative to itself
   D(i,missing(i))=0;
   fprintf('%d snp: max correlation is: %f\n', i, max(D(i,:)))
end