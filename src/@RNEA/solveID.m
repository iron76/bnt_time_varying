function obj = solveID(obj)
%solveID Inverse Dynamics with recursive Newton-Euler Algorithm (RNEA)
%   This function solves the inverse dynamics problem with the classical
%   recursive Newton-Euler algorithm, as implemented in the Featherstone's
%   toolbox. Differently from the classical implementatio, the function
%   produces an output which includes all relevant dynamic quantities. In
%   particular the output 'd' is structured as follows:
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
%   obj.IDmeas.y = [fx_1, d2q_1, ... , fx_obj.IDstate.n, d2q_obj.IDstate.n]
%
% Author: Francesco Nori
% Genova, Dec 2014

iy = 1;
for j = 1 : obj.IDsens.sensorsParams.ny
   for i = 1 : obj.IDmodel.n
      if strcmp(obj.IDsens.sensorsParams.labels{j,1}, ['fx' num2str(i)])
         obj.fx(:,i)  = obj.IDmeas.y(iy:iy+5,1);
      elseif strcmp(obj.IDsens.sensorsParams.labels{j,1}, ['d2q' num2str(i)])
         obj.d2q(i,1) = obj.IDmeas.y(iy, 1);
      end
   end
   iy = iy + obj.IDsens.sensorsParams.sizes{j,1};
end

for i = 1 : obj.IDmodel.n
   % [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
   % [~, jn{i}] = size(S{i});
   % vJ = S{i}*dq(i);
   % Xup{i} = XJ * model.Xtree{i};
   if obj.IDmodel.modelParams.parent(i) == 0
      obj.a(:,i) = obj.Xup{i}*(-obj.IDmodel.g) + obj.IDmodel.S{i}*obj.d2q(i,1);
   else
      obj.a(:,i) = obj.Xup{i}*obj.a(:,obj.IDmodel.modelParams.parent(i)) + obj.IDmodel.S{i}*obj.d2q(i,1) + crm(obj.v(:,i))*obj.vJ(:,i);
   end
   obj.fB(:,i) = obj.IDmodel.modelParams.I{i}*obj.a(:,i) + crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i);
end

obj.f = obj.fB;

if ~isempty(obj.fx) 
  for i = 1:length(obj.IDmodel.modelParams.parent)
    if ~isempty(obj.fx(:,i)) 
      obj.f(:,i) = obj.f(:,i) - obj.Xa{i}' \ obj.fx(:,i);
    end
  end
end

for i = obj.IDmodel.n:-1:1
   obj.tau(i,1) = obj.IDmodel.S{i}' * obj.f(:,i);
   if obj.IDmodel.modelParams.parent(i) ~= 0
      obj.f(:,obj.IDmodel.modelParams.parent(i)) = obj.f(:,obj.IDmodel.modelParams.parent(i)) + obj.Xup{i}'*obj.f(:,i);
   end
end

for i = 1 : obj.IDmodel.n
   obj.d((1:26)+(i-1)*26, 1) = [obj.a(:,i); obj.fB(:,i); obj.f(:,i); obj.tau(i,1); obj.fx(:,i); obj.d2q(i,1)];
end


end % solveID