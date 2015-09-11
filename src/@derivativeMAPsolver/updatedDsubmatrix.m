function [ obj ] = updatedDsubmatrix( obj )
%UPDATEDDSUBMATRIX update obj.dDdq

%% Compute the derivative of D with respect to qi (i-th element of q)
for i = 1 : obj.IDstate.n
   j = obj.IDmodel.modelParams.parent(i);
   if j > 0
      obj.dDdq{i} = set(obj.dDdq{i}, obj.dXupdq{i}, (i-1)*4+1, (j-1)*6+1);
   end
   
   for j = obj.IDmodel.sparseParams.ind_j{i}
      obj.dDdq{j} = set(obj.dDdq{j}, obj.dXupdq{j}', (i-1)*4+3, (j-1)*6+3);
   end
end

end

