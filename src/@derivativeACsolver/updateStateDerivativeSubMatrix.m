function [ obj ] = updateStateDerivativeSubMatrix( obj , d )
%UPDATESTATEDERIVATIVESUBMATRIX Compute \frac{\partial (Dd + b)}{\partial x}
%   Compute the derivative of D(q)d + b(q,\dot{q}) with respect to x
%   (x is defined as (q,\dot{q}).
%   The output is saved in the dDb (please suggest a better name) attribute
%   of the input obj.
%

for i = 1 : obj.IDmodel.n
   a{i} = d((i-1)*7+1:(i-1)*7+6, 1);
end

obj.dDb.matrix = zeros(6*obj.IDstate.n, 2*obj.IDstate.n);

%% Former non-sparse matrix
for h = 1 : obj.IDstate.n
   %% Compute D_{i,\lambda{i}}d_\lambda{i} subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      parenti = obj.IDmodel.modelParams.parent(i);
      if obj.IDmodel.modelParams.parent(i) ~= 0;
         if( h == i )
            obj.dDb = set(obj.dDb, obj.dXupdq{i} * a{parenti}, (i-1)*1+1,h);
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
         obj.dDb = set(obj.dDb, obj.dDb((i-1)*1+1,h) + ...
            obj.dXupdq{i}*(-obj.IDmodel.g) + crm(obj.dvdx{i,h})*obj.vJ(:,i), (i-1)*1+1, h);
      else 
         obj.dDb = set(obj.dDb, obj.dDb((i-1)*1+1,h) + ...
            crm(obj.dvdx{i,h})*obj.vJ(:,i), (i-1)*1+1,h);
      end
            
      % If i == 1, this term of the b vector is always zero
      if i ~= 1
         % obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
         obj.dDb = set(obj.dDb, ...
            obj.dDb((i-1)*1+1,h+obj.IDstate.n) + ...
            crm(obj.dvdx{i,obj.IDstate.n+h})*obj.vJ(:,i), (i-1)*1+1,h+obj.IDstate.n);
         % if i == h, we have also to compute the derivative of
         % obj.vJ(:,i) wrt to q_i (that is simply S_i)
         if( i == h )
            obj.dDb = set(obj.dDb, ...
               obj.dDb((i-1)*1+1,h+obj.IDstate.n) + ...
               crm(obj.v(:,i))* obj.IDmodel.S{i}, (i-1)*1+1,h+obj.IDstate.n);
         end
      end
   end
end

%% Sparse submatrix
obj.dDb_s.As = zeros(size(obj.dDb_s.As));

for h = 1 : obj.IDstate.n   
   %% Compute D_{i,\lambda{i}}d_\lambda{i} subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      parenti = obj.IDmodel.modelParams.parent(i);
      if obj.IDmodel.modelParams.parent(i) ~= 0;
         if( h == i )
            obj.dDb_s = set(obj.dDb_s, obj.dXupdq{i} * a{parenti}, (i-1)*1+1,h);
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
         obj.dDb_s = set(obj.dDb_s, obj.dDb_s((i-1)*1+1,h) + ...
            obj.dXupdq{i}*(-obj.IDmodel.g) + crm(obj.dvdx{i,h})*obj.vJ(:,i), (i-1)*1+1, h);
      else
         obj.dDb_s = set(obj.dDb_s, obj.dDb_s((i-1)*1+1,h) + ...
            crm(obj.dvdx{i,h})*obj.vJ(:,i), (i-1)*1+1,h);
      end
            
      %% db1/dvq vector
      % If i == 1, this term of the b vector is always zero
      if i ~= 1
         % obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
         obj.dDb_s = set(obj.dDb_s, ...
            obj.dDb_s((i-1)*1+1,h+obj.IDstate.n) + ...
            crm(obj.dvdx{i,obj.IDstate.n+h})*obj.vJ(:,i), (i-1)*1+1,h+obj.IDstate.n);
         % if i == h, we have also to compute the derivative of
         % obj.vJ(:,i) wrt to q_i (that is simply S_i)
         if( i == h )
            obj.dDb_s = set(obj.dDb_s, ...
               obj.dDb_s((i-1)*1+1,h+obj.IDstate.n) + ...
               crm(obj.v(:,i))* obj.IDmodel.S{i}, (i-1)*1+1,h+obj.IDstate.n);
         end
      end      
   end
end

