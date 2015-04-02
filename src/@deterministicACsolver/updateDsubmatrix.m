function [ obj ] = updateDsubmatrix( obj )
%UPDATESUBMATRIX Summary of this function goes here 	
% Detailed explanation goes here

%% Compute D_{i,i} submatrices of D matrix
% and b_{i} subvectors of b vector
for i = 1 : obj.IDstate.n
   I = (i-1)*1;
   J = (i-1)*2;
   
   if obj.IDmodel.modelParams.parent(i) == 0
      obj.b = set(obj.b, obj.Xup{i}*(-obj.IDmodel.g), I+1, 1);
   else
      obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
   end
end


%% Compute D_{i,\lambda{i}} submatrices of D matrix 
for i = 1 : obj.IDstate.n
   if obj.IDmodel.modelParams.parent(i) ~= 0
      j = obj.IDmodel.modelParams.parent(i);
      I = (i-1)*1;
      J = (j-1)*2;
      obj.D = set(obj.D, obj.Xup{i}, I+1, J+1);
   end
end


end

