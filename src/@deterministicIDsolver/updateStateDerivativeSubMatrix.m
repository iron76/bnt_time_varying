function [ obj ] = updateStateDerivativeSubMatrix( obj , d )
%UPDATESTATEDERIVATIVESUBMATRIX Compute \frac{\partial (Dd + b)}{\partial x}
%   Compute the derivative of D(q)d + b(q,\dot{q}) with respect to x
%   (x is defined as (q,\dot{q}).
%   The output is saved in the Ddbx (please suggest a better name) attribute
%   of the input obj.
%

% [a, f_B, f, tau, f_x, d2q] = extractVectors(d);

a   = ones(6, obj.IDstate.n);
f_B = ones(6, obj.IDstate.n);
f   = ones(6, obj.IDstate.n);
tau = ones(1, obj.IDstate.n);
f_x = ones(6, obj.IDstate.n);
d2q = ones(6, obj.IDstate.n);

obj.Ddbx    = zeros(19*obj.IDmodel.n, 2*obj.IDmodel.n);

for h = 1 : obj.IDstate.n
   %% Compute the derivatives of D_{i,j}d_j subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      for j = obj.IDmodel.sparseParams.ind_j{i}
         if( h == j )
            obj.Ddbx((i-1)*19+13:(i-1)*19+18,h) = (obj.dXupdq{j}')*f(:,j);
         end
      end
   end
   
   %% Compute D_{i,\lambda{i}}d_\lambda{i} subvector of Dd + b
   %  with respect to q_h (x_h)
   for i = 1 : obj.IDstate.n
      parenti = obj.IDmodel.modelParams.parent(i);
      if obj.IDmodel.modelParams.parent(i) ~= 0;
         if( h == i )
            obj.Ddbx((i-1)*19+1:(i-1)*19+6,h) = obj.dXupdq{i} * a(:,parenti);
         end
      end
   end
   
   %% Only b actually depends on \dot{q}, so we don't have to consider D
   % Compute the derivatives of D_{i,i}d_i+b_i subvector of Dd + b
   %  with respect to \dot{q}_h (x_{n+h})
   for i = 1 : obj.IDstate.n
      
      % b vector
      % If i == 1 and we have only the term of gravitational acceleration
      if (i == 1) && (h == 1)
         obj.Ddbx((i-1)*19+1:(i-1)*19+6,h) = ...
            obj.Ddbx((i-1)*19+1:(i-1)*19+6,h) + ...
            obj.dXupdq{1}*(-obj.IDmodel.g);
      elseif i ~= 1
         obj.Ddbx((i-1)*19+1:(i-1)*19+6,h) = ...
            obj.Ddbx((i-1)*19+1:(i-1)*19+6,h) + ...
            crm(obj.dvdx{i,h})*obj.vJ(:,i);
      end
      
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial q_j} =
      % \frac{\partial v_i}{\partial q_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial q_j} +
      obj.Ddbx((i-1)*19+7:(i-1)*19+12,h) = ...
         obj.Ddbx((i-1)*19+7:(i-1)*19+12,h) + ...
         crf(obj.dvdx{i,h})*obj.IDmodel.modelParams.I{i}*obj.v(:,i) + ...
         crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.dvdx{i,h};
      
      % If i == 1, this term of the b vector is always zero
      if i ~= 1
         % obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
         obj.Ddbx((i-1)*19+1:(i-1)*19+6,h+obj.IDstate.n) = ...
            obj.Ddbx((i-1)*19+1:(i-1)*19+6,h+obj.IDstate.n) + ...
            crm(obj.dvdx{i,obj.IDstate.n+h})*obj.vJ(:,i);
         % if i == h, we have also to compute the derivative of
         % obj.vJ(:,i) wrt to q_i (that is simply S_i)
         if( i == h )
            obj.Ddbx((i-1)*19+1:(i-1)*19+6,h+obj.IDstate.n) = ...
               obj.Ddbx((i-1)*19+1:(i-1)*19+6,h+obj.IDstate.n) + ...
               crm(obj.v(:,i))* obj.IDmodel.S{i};
         end
      end
      % obj.b = set(obj.b, crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i), I+2, 1);
      % \frac{\partial v_i \times^{*} I_i v_i }{\partial \dot{q}_j} =
      % \frac{\partial v_i}{\partial \dot{q}_j} \times^{*} I_i v_i +
      % v_i \times^{*} I_i \frac{\partial v_i}{\partial \dot{q}_j} +
      obj.Ddbx((i-1)*19+7:(i-1)*19+12,h+obj.IDstate.n) = ...
         obj.Ddbx((i-1)*19+7:(i-1)*19+12,h+obj.IDstate.n) + ...
         crf(obj.dvdx{i,obj.IDstate.n+h})*obj.IDmodel.modelParams.I{i}*obj.v(:,i) + ...
         crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.dvdx{i,obj.IDstate.n+h};
      
   end
end

