function [ obj ] = initdDsubmatrixIndices( obj )
%INITDDSUBMATRIX Summary of this function goes here
%   Detailed explanation goes here
iD = zeros(4*obj.IDmodel.n,1);
jD = zeros(6*obj.IDmodel.n,1);
for i = 1:obj.IDmodel.n
   iD((i-1)*4+1 : 4*i, 1) = [6 6 6 obj.IDmodel.jn(i)]';
   jD((i-1)*6+1 : 6*i, 1) = [6 6 6 obj.IDmodel.jn(i) 6 obj.IDmodel.jn(i)]';
end

I = cell(obj.IDstate.n);
J = cell(obj.IDstate.n);
for i = 1 : obj.IDstate.n
   j = obj.IDmodel.modelParams.parent(i);
   if j > 0
      % obj.dDdq{i} = set(obj.dDdq{i}, obj.dXupdq{i}, (i-1)*4+1, (j-1)*6+1);
      I{i} = [I{i}; (i-1)*4+1];
      J{i} = [J{i}; (j-1)*6+1];
   end
   
   for j = obj.IDmodel.sparseParams.ind_j{i}
      % obj.dDdq{j} = set(obj.dDdq{j}, obj.dXupdq{j}', (i-1)*4+3, (j-1)*6+3);
      I{j} = [I{j}; (i-1)*4+3];
      J{j} = [J{j}; (j-1)*6+3];
   end
end
for i = 1 : obj.IDstate.n
   obj.dDdq{i} = submatrixSparse(iD, jD, I{i}, J{i});
end   

end

