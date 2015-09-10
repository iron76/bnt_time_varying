function obj = initdDsubmatrixIndices(obj)

%% Indices for the submatrix
obj.iDb = zeros(1*obj.IDmodel.n,1);
obj.jDb = zeros(2*obj.IDmodel.n,1);
for i = 1:obj.IDmodel.n
   obj.iDb((i-1)*1+1 : 1*i, 1) = 6;
   obj.jDb((i-1)*2+1 : 2*i, 1) = [1 1]';
end

%% Indeces for the sparse matrix

obj.iDb_s = [];
obj.jDb_s = [];

for h = 1 : obj.IDstate.n
   %% Compute D_{i,\lambda{i}}d_\lambda{i} subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      % parenti = obj.IDmodel.modelParams.parent(i);
      if (obj.IDmodel.modelParams.parent(i) ~= 0) && ( h == i )         
         % obj.dDb = set(obj.dDb, obj.dXupdq{i} * a{parenti}, (i-1)*4+1,h);
         obj.iDb_s = [obj.iDb_s (i-1)*1+1];
         obj.jDb_s = [obj.jDb_s h];
      end
   end
   
   %% Only b actually depends on \dot{q}, so we don't have to consider D
   % Compute the derivatives of D_{i,i}d_i+b_i subvector of Dd + b
   %  with respect to \dot{q}_h (x_{n+h})
   for i = 1 : obj.IDstate.n
      
      % b vector
      % If i == 1 and we have only the term of gravitational acceleration
      if (i == h) && (obj.IDmodel.modelParams.parent(i) == 0)
         %obj.dDb = set(obj.dDb, obj.dDb((i-1)*4+1,h) + ...
         %   obj.dXupdq{1}*(-obj.IDmodel.g), (i-1)*4+1, h);
         obj.iDb_s = [obj.iDb_s (i-1)*1+1];
         obj.jDb_s = [obj.jDb_s h];
      elseif (i ~= h) %(i ~= h) || (obj.IDmodel.modelParams.parent(i) ~= 0)
         %obj.dDb = set(obj.dDb, obj.dDb((i-1)*4+1,h) + ...
         %   crm(obj.dvdx{i,h})*obj.vJ(:,i), (i-1)*4+1,h);
         obj.iDb_s = [obj.iDb_s (i-1)*1+1];
         obj.jDb_s = [obj.jDb_s h];
      end
      
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial q_j} =
      % \frac{\partial v_i}{\partial q_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial q_j} +
      
      % If i == 1, this term of the b vector is always zero
      if i ~= 1
         % obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
         % obj.dDb = set(obj.dDb, ...
         %   obj.dDb((i-1)*4+1,h+obj.IDstate.n) + ...
         %   crm(obj.dvdx{i,obj.IDstate.n+h})*obj.vJ(:,i), (i-1)*4+1,h+obj.IDstate.n);
         obj.iDb_s = [obj.iDb_s (i-1)*1+1];
         obj.jDb_s = [obj.jDb_s h+obj.IDstate.n];
         
         
         % if i == h, we have also to compute the derivative of
         % obj.vJ(:,i) wrt to q_i (that is simply S_i)
         if( i == h )
            % obj.dDb = set(obj.dDb, ...
            %    obj.dDb((i-1)*4+1,h+obj.IDstate.n) + ...
            %    crm(obj.v(:,i))* obj.IDmodel.S{i}, (i-1)*4+1,h+obj.IDstate.n);
            
            % Removed to avoid mutiple equivalent indexing
            % obj.iDb_s = [obj.iDb_s (i-1)*1+1];
            % obj.jDb_s = [obj.jDb_s h+obj.IDstate.n];
         end
      end
   end
end

obj.dDb_s = submatrixSparse(obj.iDb, obj.jDb, obj.iDb_s', obj.jDb_s');