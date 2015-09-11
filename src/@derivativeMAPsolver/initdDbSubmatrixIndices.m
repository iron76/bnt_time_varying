function obj = initdDbSubmatrixIndices(obj)


%% Indeces for the sparse matrix

obj.iDd_s = [];
obj.jDd_s = [];

for h = 1 : obj.IDstate.n
   %% Compute the derivatives of D_{i,j}d_j subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      for j = obj.IDmodel.sparseParams.ind_j{i}
         if( h == j )
            % obj.dDb = set(obj.dDb, (obj.dXupdq{j}')*f{j}, (i-1)*4+3,h);
            obj.iDd_s = [obj.iDd_s (i-1)*4+3];
            obj.jDd_s = [obj.jDd_s h];
         end
      end
   end
   
   %% Compute D_{i,\lambda{i}}d_\lambda{i} subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      % parenti = obj.IDmodel.modelParams.parent(i);
      if (obj.IDmodel.modelParams.parent(i) ~= 0) && ( h == i )
            % obj.dDb = set(obj.dDb, obj.dXupdq{i} * a{parenti}, (i-1)*4+1,h);
            obj.iDd_s = [obj.iDd_s (i-1)*4+1];
            obj.jDd_s = [obj.jDd_s h];
      end
   end
end
   
obj.ibD_s = [];
obj.jbD_s = [];
%% Only b actually depends on \dot{q}, so we don't have to consider D
% Compute the derivatives of D_{i,i}d_i+b_i subvector of Dd + b
%  with respect to \dot{q}_h (x_{n+h})
for h = 1 : obj.IDstate.n   
   for i = 1 : obj.IDstate.n
      
      % b vector
      % If i == 1 and we have only the term of gravitational acceleration
      if (i == h) && (obj.IDmodel.modelParams.parent(i) == 0)
         %obj.dDb = set(obj.dDb, obj.dDb((i-1)*4+1,h) + ...
         %   obj.dXupdq{1}*(-obj.IDmodel.g), (i-1)*4+1, h);
         obj.ibD_s = [obj.ibD_s (i-1)*4+1];
         obj.jbD_s = [obj.jbD_s h];
      else %(i ~= h) || (obj.IDmodel.modelParams.parent(i) ~= 0)
         %obj.dDb = set(obj.dDb, obj.dDb((i-1)*4+1,h) + ...
         %   crm(obj.dvdx{i,h})*obj.vJ(:,i), (i-1)*4+1,h);         
         obj.ibD_s = [obj.ibD_s (i-1)*4+1];
         obj.jbD_s = [obj.jbD_s h];
      end
      
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial q_j} =
      % \frac{\partial v_i}{\partial q_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial q_j} +
      
      %obj.dDb = set(obj.dDb, obj.dDb((i-1)*4+2,h) + ...
      %   crf(obj.dvdx{i,h})*obj.IDmodel.modelParams.I{i}*obj.v(:,i) + ...
      %   crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.dvdx{i,h}, (i-1)*4+2, h);
      obj.ibD_s = [obj.ibD_s (i-1)*4+2];
      obj.jbD_s = [obj.jbD_s h];      
      
      % If i == 1, this term of the b vector is always zero
      if i ~= 1
         % obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
         % obj.dDb = set(obj.dDb, ...
         %   obj.dDb((i-1)*4+1,h+obj.IDstate.n) + ...
         %   crm(obj.dvdx{i,obj.IDstate.n+h})*obj.vJ(:,i), (i-1)*4+1,h+obj.IDstate.n);
         obj.ibD_s = [obj.ibD_s (i-1)*4+1];
         obj.jbD_s = [obj.jbD_s h+obj.IDstate.n];
         
         
         % if i == h, we have also to compute the derivative of
         % obj.vJ(:,i) wrt to q_i (that is simply S_i)
         if( i == h ) 
            % obj.dDb = set(obj.dDb, ...
            %    obj.dDb((i-1)*4+1,h+obj.IDstate.n) + ...
            %    crm(obj.v(:,i))* obj.IDmodel.S{i}, (i-1)*4+1,h+obj.IDstate.n);
            
            % Removed to avoid mutiple equivalent indexing
            % obj.ibD_s = [obj.ibD_s (i-1)*4+1];
            % obj.jbD_s = [obj.jbD_s h+obj.IDstate.n];
         end
      end
      % obj.b = set(obj.b, crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i), I+2, 1);
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial \dot{q}_j} =
      % \frac{\partial v_i}{\partial \dot{q}_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial \dot{q}_j} +
      % obj.dDb = set(obj.dDb, ...
      %    obj.dDb((i-1)*4+2,h+obj.IDstate.n) + ...
      %    crf(obj.dvdx{i,obj.IDstate.n+h})*obj.IDmodel.modelParams.I{i}*obj.v(:,i) + ...
      %    crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.dvdx{i,obj.IDstate.n+h}, (i-1)*4+2,h+obj.IDstate.n);
      obj.ibD_s = [obj.ibD_s (i-1)*4+2];
      obj.jbD_s = [obj.jbD_s h+obj.IDstate.n];
   end
end

ibD = zeros(4*obj.IDmodel.n,1);
jbD = zeros(2*obj.IDmodel.n,1);
for i = 1:obj.IDmodel.n
   ibD((i-1)*4+1 : 4*i, 1) = [6 6 6 obj.IDmodel.jn(i)]';
   jbD((i-1)*2+1 : 2*i, 1) = [1 1];
end
obj.dbD_s = submatrixSparse(ibD, jbD, obj.ibD_s', obj.jbD_s');
obj.dDd_s = submatrixSparse(ibD, jbD, obj.iDd_s', obj.jDd_s');
