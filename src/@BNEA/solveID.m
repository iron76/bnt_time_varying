function obj = solveID(obj)
%solveID Inverse Dynamics with Bayesian networks (BNEA)
%   This function solves the inverse dynamics problem with the a Bayesian
%   network representation, as described in the paper "BERDY: Bayesian
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
%   matrix Y(q, dq), is represented as a Bayesian network, implemented
%   using the Matlab toolbox provided by Kevin Muprhy. Moreover, the
%   variables d should satisfy the Newton-Euler equations represented as
%   D(q,dq) d + b(q, dq) = 0, again represented as a Bayesian network.
%
% Author: Francesco Nori
% Genova, Dec 2014

NB     = obj.IDmodel.modelParams.NB;

% % for i = 1:obj.IDmodel.modelParams.NB
% %   [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
% %   [~, jn{i}] = size(S{i});
% %   for j = 1:obj.IDmodel.modelParams.NB
% %     Dc{i,j} = zeros(18+jn{i}, 24+2*jn{i});
% %   end
% %   vJ = S{i}*qd(i);
% %   Xup{i} = XJ * model.Xtree{i};
% %   if model.parent(i) == 0
% %     v{i} = vJ;
% %     % a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
% %   else
% %     v{i} = Xup{i}*v{model.parent(i)} + vJ;
% %     % a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
% %   end
% % end
% %
% % if nargin >= 5
% %   [~, Xa] = apply_external_forces( model.parent, Xup, cell(NB), cell(NB) );
% % end


for i = 1:obj.IDmodel.modelParams.NB
   if obj.IDmodel.modelParams.parent(i) == 0
      % a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
      ai   = obj.bnt.nodes.index{i}(1);
      nai  = obj.bnt.nodes.sizes{i,1}(1);
      d2qi = obj.bnt.nodes.index{i}(6);
      pars  = cell2mat(obj.bnt.bnet.parents(ai));
      Wa    = cell(1, length(pars));
      for j = 1 : length(pars)
         if pars(j) == d2qi
            Wa{1, j} = obj.IDmodel.S{i};
         end
      end
      W = cell2mat(Wa);
      obj.bnt.bnet.CPD{ai} = gaussian_CPD(obj.bnt.bnet, ai, 'mean', obj.Xup{i}*(-obj.IDmodel.g), 'cov', obj.IDmodel.modelParams.Sm.a{i}, 'weights', W);
   else
      % a{i} = ... + S{i}*qdd(i) + crm(v{i})*vJ;
      % a{i} = Xup{i}*a{model.parent(i)} + ...
      ai      = obj.bnt.nodes.index{i}(1);
      d2qi    = obj.bnt.nodes.index{i}(6);
      nai     = obj.bnt.nodes.sizes{i,1}(1);
      aj      = obj.bnt.nodes.index{obj.IDmodel.modelParams.parent(i)}(1);
      
      pars  = cell2mat(obj.bnt.bnet.parents(ai));
      Wa    = cell(1, length(pars));
      for j = 1 : length(pars)
         if pars(j) == aj
            Wa{1, j} = obj.Xup{i};
         elseif pars(j) == d2qi
            Wa{1, j} = obj.IDmodel.S{i};
         end
      end
      W = cell2mat(Wa);
      obj.bnt.bnet.CPD{ai} = gaussian_CPD(obj.bnt.bnet, ai, 'mean', crm(obj.v(:,i))*obj.vJ(:,i), 'cov', obj.IDmodel.modelParams.Sm.a{i}, 'weights', W);
   end
   % fB{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
   fBi  = obj.bnt.nodes.index{i}(2);
   nfBi = obj.bnt.nodes.sizes{i,1}(2);
   pars  = cell2mat(obj.bnt.bnet.parents(fBi));
   Wa    = cell(1, length(pars));
   for j = 1 : length(pars)
      if pars(j) == ai
         Wa{1, j} = obj.IDmodel.modelParams.I{i};
      end
   end
   W = cell2mat(Wa);
   obj.bnt.bnet.CPD{fBi} = gaussian_CPD(obj.bnt.bnet, fBi, 'mean', crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i), 'cov', obj.IDmodel.modelParams.Sm.fB{i}, 'weights', W);
end

for i = obj.IDmodel.modelParams.NB:-1:1
   % f{i} = fB{i} - Xa{i}' \ f_ext{i};
   % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};
   fBi = obj.bnt.nodes.index{i}(2);
   fi  = obj.bnt.nodes.index{i}(3);
   fxi = obj.bnt.nodes.index{i}(5);
   nfi = obj.bnt.nodes.sizes{i,1}(3);
   pars  = cell2mat(obj.bnt.bnet.parents(fi));
   Wa    = cell(1, length(pars));
   for j = 1 : length(pars)
      if pars(j) == fBi
         Wa{1, j} = eye(nfi);
      elseif pars(j) == fxi
         Wa{1, j} = -inv(obj.Xa{i}');
      else
         Wa{1, j} = obj.Xup{obj.bnt.link(pars(j))}';
      end
   end
   W = cell2mat(Wa);
   
   obj.bnt.bnet.CPD{fi} = gaussian_CPD(obj.bnt.bnet, fi, 'mean', zeros(nfi,1), 'cov', obj.IDmodel.modelParams.Sm.f{i}, 'weights', W);
   
   % tau(i,1) = S{i}' * f{i};
   taui  = obj.bnt.nodes.index{i}(4);
   ntaui = obj.bnt.nodes.sizes{i,1}(4);
   pars  = cell2mat(obj.bnt.bnet.parents(taui));
   Wa    = cell(1, length(pars));
   for j = 1 : length(pars)
      if pars(j) == fi
         Wa{1, j} = obj.IDmodel.S{i}';
      end
   end
   W = cell2mat(Wa);
   obj.bnt.bnet.CPD{taui} = gaussian_CPD(obj.bnt.bnet, taui, 'mean', zeros(ntaui, 1), 'cov', obj.IDmodel.modelParams.Sm.tau{i}, 'weights', W);
   
   % fxi
   fxi  = obj.bnt.nodes.index{i}(5);
   nfxi = obj.bnt.nodes.sizes{i,1}(5);
   obj.bnt.bnet.CPD{fxi} = gaussian_CPD(obj.bnt.bnet, fxi, 'mean', zeros(nfxi, 1), 'cov', obj.IDmodel.modelParams.Su.fx{i});
   
   % d2qi
   d2qi  = obj.bnt.nodes.index{i}(6);
   nd2qi = obj.bnt.nodes.sizes{i,1}(6);
   obj.bnt.bnet.CPD{d2qi} = gaussian_CPD(obj.bnt.bnet, d2qi, 'mean', zeros(nd2qi, 1), 'cov', obj.IDmodel.modelParams.Su.d2q{i});
   
end

for i = 1 : obj.IDsens.sensorsParams.ny
   yi  = obj.bnt.nodes.index{NB*6 + i};
   nyi = obj.bnt.nodes.sizes{NB*6 + i};
   obj.bnt.bnet.CPD{yi} = gaussian_CPD(obj.bnt.bnet, yi, 'mean', zeros(nyi, 1), 'cov', obj.IDsens.sensorsParams.Sy{i,1}, 'weights', obj.bnt.Wy{i,1});
end

a   = zeros(NB,6);
fB  = zeros(NB,6);
f   = zeros(NB,6);
tau = zeros(NB,1);

Sa   = cell(NB,1);
SfB  = cell(NB,1);
Sf   = cell(NB,1);
Stau = cell(NB,1);
Sfx  = cell(NB,1);
Sd2q = cell(NB,1);

if isempty(obj.bnt.engine)
   obj.bnt.engine = jtree_inf_engine(obj.bnt.bnet);
else
   % disp('Using previous path on jtree!')
   obj.bnt.engine = bnet_to_engine(obj.bnt.bnet, obj.bnt.engine);
end


evidence  = cell(1,6*NB+obj.IDsens.sensorsParams.ny);
iy = 1;
for i = 1 : obj.IDsens.sensorsParams.ny
   evidence{obj.bnt.nodes.index{NB*6 + i}} = obj.IDmeas.y(iy:(iy+obj.IDsens.sensorsParams.sizes{i,1}-1));
   iy = iy + obj.IDsens.sensorsParams.sizes{i,1};
end
obj.bnt.engine = enter_evidence(obj.bnt.engine, evidence);
for i = NB:-1:1
   tmp        = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(1));
   a(i,:)   = tmp.mu';
   Sa{i}    = tmp.Sigma;
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(2));
   fB(i,:)  = tmp.mu';
   SfB{i}   = tmp.Sigma;
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(3));
   f(i,:)   = tmp.mu';
   Sf{i}    = tmp.Sigma;
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(4));
   tau(i,1) = tmp.mu;
   Stau{i}  = tmp.Sigma;
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(5));
   fx(i,:)  = tmp.mu';
   Sfx{i}   = tmp.Sigma;
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(6));
   d2q(i,1) = tmp.mu';
   Sd2q{i}  = tmp.Sigma;
   
end

for i = 1 : NB
   obj.Sd((1:6)+(i-1)*19, (1:6)+(i-1)*19) = Sa{i};
   obj.Sd((7:12)+(i-1)*19, (7:12)+(i-1)*19) = SfB{i};
   obj.Sd((13:18)+(i-1)*19, (13:18)+(i-1)*19) = Sf{i};
   obj.Sd((19:19)+(i-1)*19, (19:19)+(i-1)*19) = Stau{i};
   
   obj.d((1:26)+(i-1)*26, 1) = [a( i,1:6), fB(i,1:6), f( i,1:6), tau(i,1), fx(i,1:6), d2q(i,1)]';
end

for i = 1 : NB
   obj.Sd((1:6)+(i-1)*7+19*NB, (1:6)+(i-1)*7+19*NB) = Sfx{i};
   obj.Sd((7:7)+(i-1)*7+19*NB, (7:7)+(i-1)*7+19*NB) = Sd2q{i};
   
   obj.d((1:26)+(i-1)*26, 1) = [a( i,1:6), fB(i,1:6), f( i,1:6), tau(i,1), fx(i,1:6), d2q(i,1)]';
end

end % solveID