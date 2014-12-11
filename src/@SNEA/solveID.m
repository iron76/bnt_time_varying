function obj = solveID(obj)
%solveID Inverse Dynamics with sparse Newton-Euler Algorithm (SNEA)
%   This function solves the inverse dynamics problem with the sparse
%   Newton-Euler algorithm, as described in the paper "BERDY: Bayesian 
%   Estimation for Robot Dynamics. A Probabilistic Estimation of Whole-Body
%   Dynamics with Redundant Measurements." The output 'd' is structured as
%   follows:
%
%   d   = [d_1, d_2, ..., d_obj.IDstate.n]
%
%   where:
%
%   d_i = [a_i, fB_i, f_i, tau_i, fx_i, d2q_i]
%
%   and a_i is the link-i spatial accelration, fB_i is the net spatial
%   force on the link-i, f_i is spatial wrench transmitted to link-i from
%   its parent, tau_i is torque on joint-i, fx_i is the external force on
%   link-i and d2q_i is acceleration of joint-i. The input to the algorithm
%   is in obj.IDmeas.y organized as follows:
%
%   obj.IDmeas.y = [y_1, y_2, ... , y_obj.IDsens.m]
%
%   The relationship between d and y is given by Y(q, dq) d = y where the
%   matrix Y(q, dq), is represented as a sparse matrix. Moreover, the
%   variables d should satisfy the Newton-Euler equations represented as
%   D(q,dq) d + b(q, dq) = 0, again represented as a sparse matrix. 
%
% Author: Francesco Nori
% Genova, Dec 2014

for i = 1 : obj.IDmodel.n
   obj.fx(:,i)  = obj.IDmeas.y((1:6)+(i-1)*7,1);
   obj.d2q(i,1) = obj.IDmeas.y(7*i);
end

a_grav = get_gravity(obj.IDmodel.modelParams);

for i = 1 : obj.IDmodel.n
   % [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
   % [~, jn{i}] = size(S{i});
   % vJ = S{i}*dq(i);
   % Xup{i} = XJ * model.Xtree{i};
   if obj.IDmodel.modelParams.parent(i) == 0
      obj.v(:,i) = obj.vJ(:,i);
      obj.a(:,i) = obj.Xup{i}*(-a_grav) + obj.IDmodel.S(:,i)*obj.d2q(i,1);
   else
      obj.v(:,i) = obj.Xup{i}*obj.v(:,obj.IDmodel.modelParams.parent(i)) + obj.vJ(:,i);
      obj.a(:,i) = obj.Xup{i}*obj.a(:,obj.IDmodel.modelParams.parent(i)) + obj.IDmodel.S(:,i)*obj.d2q(i,1) + crm(obj.v(:,i))*obj.vJ(:,i);
   end
   obj.fB(:,i) = obj.IDmodel.modelParams.I{i}*obj.a(:,i) + crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i);
end

obj.f = obj.fB;

if ~isempty(obj.fx) 
  for i = 1:length(obj.IDmodel.modelParams.parent)
    if obj.IDmodel.modelParams.parent(i) == 0
      obj.Xa{i} = obj.Xup{i};
    else
      obj.Xa{i} = obj.Xup{i} * obj.Xa{obj.IDmodel.modelParams.parent(i)};
    end
    if ~isempty(obj.fx(:,i)) 
      obj.f(:,i) = obj.f(:,i) - obj.Xa{i}' \ obj.fx(:,i);
    end
  end
end

for i = obj.IDmodel.n:-1:1
   obj.tau(i,1) = obj.IDmodel.S(:,i)' * obj.f(:,i);
   if obj.IDmodel.modelParams.parent(i) ~= 0
      obj.f(:,obj.IDmodel.modelParams.parent(i)) = obj.f(:,obj.IDmodel.modelParams.parent(i)) + obj.Xup{i}'*obj.f(:,i);
   end
end

for i = 1 : obj.IDmodel.n
   obj.d((1:26)+(i-1)*26, 1) = [obj.a(:,i); obj.fB(:,i); obj.f(:,i); obj.tau(i,1); obj.fx(:,i); obj.d2q(i,1)];
end


end % solveID