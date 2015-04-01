function [ obj ] = updateSubMatrix( obj )
%UPDATESUBMATRIX Summary of this function goes here
%   Detailed explanation goes here
for i = 1 : obj.IDstate.n
   I = (i-1)*4;
   J = (i-1)*6;
   
   if i == 1
      obj.b = set(obj.b, obj.Xup{i}*(-obj.IDmodel.g), I+1, 1);
   else
      obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
   end
   obj.b = set(obj.b, crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i), I+2, 1);
   obj.D = set(obj.D, -inv(obj.Xa{i}'), I+3, J+5);
end

for i = 1 : obj.IDstate.n
   for j = obj.IDmodel.sparseParams.ind_j{i}
      I = (i-1)*4;
      J = (j-1)*6;
      % f{obj.IDmodel.modelParams.parent(j)} = f{obj.IDmodel.modelParams.parent(j)} + obj.Xup{j}'*f{j};
      % Dc{i,j} = [ zeros(12, 24+2*obj.jn(i))
      %     zeros(6,6) zeros(6,6) obj.Xup{j}' zeros(6, obj.jn(i)) zeros(6,6) zeros(6, obj.jn(i))
      %     zeros(obj.jn(i), 24+2*obj.jn(i))];
      obj.D = set(obj.D, obj.Xup{j}', I+3, J+3);
   end
end

for i = 1 : obj.IDstate.n
   if obj.IDmodel.modelParams.parent(i) ~= 0
      j = obj.IDmodel.modelParams.parent(i);
      I = (i-1)*4;
      J = (j-1)*6;
      obj.D = set(obj.D, obj.Xup{i}, I+1, J+1);
   end
end


end

