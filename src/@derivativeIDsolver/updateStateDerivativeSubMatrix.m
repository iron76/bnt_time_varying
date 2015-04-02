function [ obj ] = updateStateDerivativeSubMatrix( obj , d )
%UPDATESTATEDERIVATIVESUBMATRIX Compute \frac{\partial (Dd + b)}{\partial x}
%   Compute the derivative of D(q)d + b(q,\dot{q}) with respect to x
%   (x is defined as (q,\dot{q}).
%   The output is saved in the dDb (please suggest a better name) attribute
%   of the input obj.
%

[a, ~, f, ~, ~, ~] = extractDynVar(obj.IDmodel.n, d);

obj.dDb.matrix = zeros(19*obj.IDstate.n, 2*obj.IDstate.n);

%% Former non-sparse matrix
for h = 1 : obj.IDstate.n
   %% Compute the derivatives of D_{i,j}d_j subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      for j = obj.IDmodel.sparseParams.ind_j{i}
         if( h == j )
            obj.dDb = set(obj.dDb, (obj.dXupdq{j}')*f{j}, (i-1)*4+3,h);
         end
      end
   end
   
   %% Compute D_{i,\lambda{i}}d_\lambda{i} subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      parenti = obj.IDmodel.modelParams.parent(i);
      if obj.IDmodel.modelParams.parent(i) ~= 0;
         if( h == i )
            obj.dDb = set(obj.dDb, obj.dXupdq{i} * a{parenti}, (i-1)*4+1,h);
         end
      end
   end
   
   %% Only b actually depends on \dot{q}, so we don't have to consider D
   % Compute the derivatives of D_{i,i}d_i+b_i subvector of Dd + b
   %  with respect to \dot{q}_h (x_{n+h})
   for i = 1 : obj.IDstate.n
      
      % b vector
      % If parent(i) == 1 and we have only the term of gravitational acceleration
      if (obj.IDmodel.modelParams.parent(i)==0) && (h == i)
         obj.dDb = set(obj.dDb, obj.dDb((i-1)*4+1,h) + ...
            obj.dXupdq{i}*(-obj.IDmodel.g) + crm(obj.dvdx{i,h})*obj.vJ(:,i), (i-1)*4+1, h);
      else 
         obj.dDb = set(obj.dDb, obj.dDb((i-1)*4+1,h) + ...
            crm(obj.dvdx{i,h})*obj.vJ(:,i), (i-1)*4+1,h);
      end
      
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial q_j} =
      % \frac{\partial v_i}{\partial q_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial q_j} +
      obj.dDb = set(obj.dDb, obj.dDb((i-1)*4+2,h) + ...
         crf(obj.dvdx{i,h})*obj.IDmodel.modelParams.I{i}*obj.v(:,i) + ...
         crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.dvdx{i,h}, (i-1)*4+2, h);
      
      % If i == 1, this term of the b vector is always zero
      if i ~= 1
         % obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
         obj.dDb = set(obj.dDb, ...
            obj.dDb((i-1)*4+1,h+obj.IDstate.n) + ...
            crm(obj.dvdx{i,obj.IDstate.n+h})*obj.vJ(:,i), (i-1)*4+1,h+obj.IDstate.n);
         % if i == h, we have also to compute the derivative of
         % obj.vJ(:,i) wrt to q_i (that is simply S_i)
         if( i == h )
            obj.dDb = set(obj.dDb, ...
               obj.dDb((i-1)*4+1,h+obj.IDstate.n) + ...
               crm(obj.v(:,i))* obj.IDmodel.S{i}, (i-1)*4+1,h+obj.IDstate.n);
         end
      end
      % obj.b = set(obj.b, crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i), I+2, 1);
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial \dot{q}_j} =
      % \frac{\partial v_i}{\partial \dot{q}_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial \dot{q}_j} +
      obj.dDb = set(obj.dDb, ...
         obj.dDb((i-1)*4+2,h+obj.IDstate.n) + ...
         crf(obj.dvdx{i,obj.IDstate.n+h})*obj.IDmodel.modelParams.I{i}*obj.v(:,i) + ...
         crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.dvdx{i,obj.IDstate.n+h}, (i-1)*4+2,h+obj.IDstate.n);
      
   end
end

%% Sparse submatrix
obj.dDb_s.As = zeros(size(obj.dDb_s.As));

for h = 1 : obj.IDstate.n
   %% Compute the derivatives of D_{i,j}d_j subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      for j = obj.IDmodel.sparseParams.ind_j{i}
         if( h == j )
            obj.dDb_s = set(obj.dDb_s, (obj.dXupdq{j}')*f{j}, (i-1)*4+3,h);
         end
      end
   end
   
   %% Compute D_{i,\lambda{i}}d_\lambda{i} subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      parenti = obj.IDmodel.modelParams.parent(i);
      if obj.IDmodel.modelParams.parent(i) ~= 0;
         if( h == i )
            obj.dDb_s = set(obj.dDb_s, obj.dXupdq{i} * a{parenti}, (i-1)*4+1,h);
         end
      end
   end
   
   %% Only b actually depends on \dot{q}, so we don't have to consider D
   % Compute the derivatives of D_{i,i}d_i+b_i subvector of Dd + b
   %  with respect to \dot{q}_h (x_{n+h})
   for i = 1 : obj.IDstate.n
      
      
      
      %% db1/dq vector
      % If parent(i) == 0 and we have only the term of gravitational acceleration
      if (obj.IDmodel.modelParams.parent(i)==0) && (h == i)
         obj.dDb_s = set(obj.dDb_s, obj.dDb_s((i-1)*4+1,h) + ...
            obj.dXupdq{i}*(-obj.IDmodel.g) + crm(obj.dvdx{i,h})*obj.vJ(:,i), (i-1)*4+1, h);
      else
         obj.dDb_s = set(obj.dDb_s, obj.dDb_s((i-1)*4+1,h) + ...
            crm(obj.dvdx{i,h})*obj.vJ(:,i), (i-1)*4+1,h);
      end
      
      %% db2/dq vector
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial q_j} =
      % \frac{\partial v_i}{\partial q_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial q_j} +
      obj.dDb_s = set(obj.dDb_s, obj.dDb_s((i-1)*4+2,h) + ...
         crf(obj.dvdx{i,h})*obj.IDmodel.modelParams.I{i}*obj.v(:,i) + ...
         crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.dvdx{i,h}, (i-1)*4+2, h);
      
      %% db1/dvq vector
      % If i == 1, this term of the b vector is always zero
      if i ~= 1
         % obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
         obj.dDb_s = set(obj.dDb_s, ...
            obj.dDb_s((i-1)*4+1,h+obj.IDstate.n) + ...
            crm(obj.dvdx{i,obj.IDstate.n+h})*obj.vJ(:,i), (i-1)*4+1,h+obj.IDstate.n);
         % if i == h, we have also to compute the derivative of
         % obj.vJ(:,i) wrt to q_i (that is simply S_i)
         if( i == h )
            obj.dDb_s = set(obj.dDb_s, ...
               obj.dDb_s((i-1)*4+1,h+obj.IDstate.n) + ...
               crm(obj.v(:,i))* obj.IDmodel.S{i}, (i-1)*4+1,h+obj.IDstate.n);
         end
      end
      
      %% db2/dvq vector
      % obj.b = set(obj.b, crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i), I+2, 1);
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial \dot{q}_j} =
      % \frac{\partial v_i}{\partial \dot{q}_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial \dot{q}_j} +
      obj.dDb_s = set(obj.dDb_s, ...
         obj.dDb_s((i-1)*4+2,h+obj.IDstate.n) + ...
         crf(obj.dvdx{i,obj.IDstate.n+h})*obj.IDmodel.modelParams.I{i}*obj.v(:,i) + ...
         crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.dvdx{i,obj.IDstate.n+h}, (i-1)*4+2,h+obj.IDstate.n);
      
   end
end

