function [ obj ] = updatedDbSubmatrix( obj , d )
%UPDATESTATEDERIVATIVESUBMATRIX Compute \frac{\partial (Dd + b)}{\partial x}
%   Compute the derivative of D(q)d + b(q,\dot{q}) with respect to x
%   (x is defined as (q,\dot{q}).
%   The output is saved in the dDb (please suggest a better name) attribute
%   of the input obj.
%

[m,n] = size(d);
if (m ~= obj.IDmodel.modelParams.NB * 26) || (n ~= 1)
   error('[ERROR] The input d should be provided as a column vector with 26*model.NB rows');
end

[a, ~, f, ~, ~, ~] = extractDynVar(obj.IDmodel.n, d);

%% Sparse submatrix
obj.dDd_s.As = zeros(size(obj.dDd_s.As));

ks = obj.dDd_s.ks;
ps = obj.dDd_s.ps;
As = obj.dDd_s.As;

for h = 1 : obj.IDstate.n
   %% Compute the derivatives of D_{i,j}d_j subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      for j = obj.IDmodel.sparseParams.ind_j{i}
         if( h == j )
            k = ks((i-1)*4+3,h);
            A = (obj.dXupdq{j}')*f{j};
            % obj.dbD_s.As(ps(k)+1:ps(k+1)) = A(:);
            As(ps(k)+1:ps(k+1)) = A(:);
            
            %obj.dDb = set(obj.dDb, (obj.dXupdq{j}')*f{j}, (i-1)*4+3,h);
         end
      end
   end
   
   %% Compute D_{i,\lambda{i}}d_\lambda{i} subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      parenti = obj.IDmodel.modelParams.parent(i);
      if obj.IDmodel.modelParams.parent(i) ~= 0;
         if( h == i )
            k = ks((i-1)*4+1,h);
            A = obj.dXupdq{i} * a{parenti};
            % obj.dbD_s.As(ps(k)+1:ps(k+1)) = A(:);
            As(ps(k)+1:ps(k+1)) = A(:);
            
            % obj.dbD_s = set(obj.dbD_s, obj.dXupdq{i} * a{parenti}, (i-1)*4+1,h);
         end
      end
   end
end
obj.dDd_s.As = As;


obj.dbD_s.As = zeros(size(obj.dbD_s.As));

ks = obj.dbD_s.ks;
ps = obj.dbD_s.ps;

is = obj.dbD_s.is;
js = obj.dbD_s.js;
As = obj.dbD_s.As;

for h = 1 : obj.IDstate.n
   %% Only b actually depends on \dot{q}, so we don't have to consider D
   % Compute the derivatives of D_{i,i}d_i+b_i subvector of Dd + b
   %  with respect to \dot{q}_h (x_{n+h})
   for i = 1 : obj.IDstate.n
      %% db1/dq vector
      % If parent(i) == 0 and we have only the term of gravitational acceleration
      if (obj.IDmodel.modelParams.parent(i)==0) && (h == i)
         k = ks((i-1)*4+1, h);
         I = is(ps(k)+1);
         J = js(ps(k)+1);
         l = ps(k)+1:ps(k+1);
         B = sparse(is(l) - I + 1, js(l) - J + 1, As(l));
         
         A = B + obj.dXupdq{i}*(-obj.IDmodel.g);
         % obj.dbD_s.As(ps(k)+1:ps(k+1)) = A(:);
         As(ps(k)+1:ps(k+1)) = A(:);
         
         %obj.dbD_s = set(obj.dbD_s, obj.dbD_s((i-1)*4+1,h) + ...
         %   obj.dXupdq{i}*(-obj.IDmodel.g) + crmmult(obj.dvdx{i,h},obj.vJ(:,i)), (i-1)*4+1, h);
      else
         k = ks((i-1)*4+1,h);
         I = is(ps(k)+1);
         J = js(ps(k)+1);
         l = ps(k)+1:ps(k+1);
         B = sparse(is(l) - I + 1, js(l) - J + 1, As(l));
         
         A = B + crmmult(obj.dvdx{i,h}, obj.vJ(:,i));
         % obj.dbD_s.As(ps(k)+1:ps(k+1)) = A(:);
         As(ps(k)+1:ps(k+1)) = A(:);
         % obj.dbD_s = set(obj.dbD_s, obj.dbD_s((i-1)*4+1,h) + ...
         %   crmmult(obj.dvdx{i,h}, obj.vJ(:,i)), (i-1)*4+1,h);
      end
      
      %% db2/dq vector
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial q_j} =
      % \frac{\partial v_i}{\partial q_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial q_j} +
      k = ks((i-1)*4+2, h);
      I = is(ps(k)+1);
      J = js(ps(k)+1);
      l = ps(k)+1:ps(k+1);
      B = sparse(is(l) - I + 1, js(l) - J + 1, As(l));
      
      A = B + crfmult(obj.dvdx{i,h}, obj.IDmodel.modelParams.I{i}*obj.v(:,i)) + crfmult(obj.v(:,i), obj.IDmodel.modelParams.I{i}*obj.dvdx{i,h});
      % obj.dbD_s.As(ps(k)+1:ps(k+1)) = A(:);
      As(ps(k)+1:ps(k+1)) = A(:);
      
      
      %obj.dbD_s = set(obj.dbD_s, obj.dbD_s((i-1)*4+2,h) + ...
      %   crfmult(obj.dvdx{i,h}, obj.IDmodel.modelParams.I{i}*obj.v(:,i)) + ...
      %   crfmult(obj.v(:,i), obj.IDmodel.modelParams.I{i}*obj.dvdx{i,h}), (i-1)*4+2, h);
      
      %% db1/dvq vector
      % If i == 1, this term of the b vector is always zero
      if i ~= 1
         % obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
         k = ks((i-1)*4+1,h+obj.IDstate.n);
         I = is(ps(k)+1);
         J = js(ps(k)+1);
         l = ps(k)+1:ps(k+1);
         B = sparse(is(l) - I + 1, js(l) - J + 1, As(l));
         
         A = B + crmmult(obj.dvdx{i,obj.IDstate.n+h}, obj.vJ(:,i));
         % obj.dbD_s.As(ps(k)+1:ps(k+1)) = A(:);
         As(ps(k)+1:ps(k+1)) = A(:);
         
         %obj.dbD_s = set(obj.dbD_s, ...
         %   obj.dbD_s((i-1)*4+1,h+obj.IDstate.n) + ...
         %   crmmult(obj.dvdx{i,obj.IDstate.n+h}, obj.vJ(:,i)), (i-1)*4+1,h+obj.IDstate.n);
         % if i == h, we have also to compute the derivative of
         % obj.vJ(:,i) wrt to q_i (that is simply S_i)
         if( i == h )
            k = ks((i-1)*4+1,h+obj.IDstate.n);
            I = is(ps(k)+1);
            J = js(ps(k)+1);
            l = ps(k)+1:ps(k+1);
            B = sparse(is(l) - I + 1, js(l) - J + 1, As(l));
            
            A = B + crmmult(obj.v(:,i), obj.IDmodel.S{i});
            % obj.dbD_s.As(ps(k)+1:ps(k+1)) = A(:);
            As(ps(k)+1:ps(k+1)) = A(:);
            
            %obj.dbD_s = set(obj.dbD_s, ...
            %   obj.dbD_s((i-1)*4+1,h+obj.IDstate.n) + ...
            %   crmmult(obj.v(:,i), obj.IDmodel.S{i}), (i-1)*4+1,h+obj.IDstate.n);
         end
      end
      
      %% db2/dvq vector
      % obj.b = set(obj.b, crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i), I+2, 1);
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial \dot{q}_j} =
      % \frac{\partial v_i}{\partial \dot{q}_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial \dot{q}_j} +
      k = ks((i-1)*4+2,h+obj.IDstate.n);
      I = is(ps(k)+1);
      J = js(ps(k)+1);
      l = ps(k)+1:ps(k+1);
      B = sparse(is(l) - I + 1, js(l) - J + 1, As(l));
      
      A = B + crfmult(obj.dvdx{i,obj.IDstate.n+h}, obj.IDmodel.modelParams.I{i}*obj.v(:,i)) + crfmult(obj.v(:,i), obj.IDmodel.modelParams.I{i}*obj.dvdx{i,obj.IDstate.n+h});
      % obj.dbD_s.As(ps(k)+1:ps(k+1)) = A(:);
      As(ps(k)+1:ps(k+1)) = A(:);
      
      % obj.dbD_s = set(obj.dbD_s, ...
      %   obj.dbD_s((i-1)*4+2,h+obj.IDstate.n) + ...
      %   crfmult(obj.dvdx{i,obj.IDstate.n+h}, obj.IDmodel.modelParams.I{i}*obj.v(:,i)) + ...
      %   crfmult(obj.v(:,i), obj.IDmodel.modelParams.I{i}*obj.dvdx{i,obj.IDstate.n+h}), (i-1)*4+2,h+obj.IDstate.n);
      
   end
end
obj.dbD_s.As = As;
