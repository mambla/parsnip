numOfSNP=300;
thresh = 0.8;

dist = zero(i,3);
for i=1:numOfSNP
   tmp = D(i,:);
   a = sort(tmp, 'descend');
   a = a(1:3);
   
   %ignore snp with top corelated < thresh
   if a(1) < thresh
       continue;
   end
   
   %the real index of this snp
   ind = find(tmp==0);
   for j=1:3
       %find the index of snp that is the (1/2/3) most correlated with #i
       indexes = find(tmp==a(j));       
       %find is distance from #i
       dist(i,j) = min(abs(repmat([ind],size(indexes))-indexes));
   end
   
   fprintf('snp #%d:\t%.2f@%d\t%.2f@%d\t%.2f@%d\n',...
       i,...
       a(1), dist(i,1),...
       a(2), dist(i,2),...
       a(3), dist(i,3));
end