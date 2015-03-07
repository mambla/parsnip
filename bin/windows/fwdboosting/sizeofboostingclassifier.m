function s = sizeofboostingclassifier(carcl,NUMBOOSTINGROUNDS);

NUMBOOSTINGROUNDS = min(NUMBOOSTINGROUNDS,length(carcl));
tmp = zeros(NUMBOOSTINGROUNDS,1);
for i = 1:NUMBOOSTINGROUNDS; tmp(i) = carcl(i).featureNdx; end
featruestocompute = unique(tmp);
s = length(featruestocompute);
