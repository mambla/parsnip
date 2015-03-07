function [labels,weights] = CLSgentleBoostC(X,Model);

if ~Model.Nbagging 
  [Cx,Fx] = strongGentleClassifier(X,Model.classifier); 
  weights = Fx;
  Cx = (Cx+3)/2;
else
  Cx = [];
  Fx = [];
  for i = 1:Model.Nbagging,
    [tCx,tFx] = strongGentleClassifier(X,Model.classifier{i});
    tCx = tCx';
    tFx = tFx';
    Cx = [Cx, tCx];
    Fx = [Fx, tFx];
  end
 
  if 0 
    weights = mean(Fx')';
  else
    weights = mean(Cx')';
  end
  Cx = sign(sum(Cx')');
  Cx = (Cx+3)/2;
end
labels = Model.uniquey(floor(Cx));

