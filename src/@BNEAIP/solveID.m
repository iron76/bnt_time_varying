function obj = solveID(obj)
%solveID Inverse Dynamics with Bayesian networks with Inertial Parameters (BNEAIP)
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
%   using the Matlab toolbox provided by Kevin Murohy. Moreover, the
%   variables d should satisfy the Newton-Euler equations represented as
%   D(q,dq) d + b(q, dq) = 0, again represented as a Bayesian network.
%
%   BNEAIP is like the BNEA, but with additional support for inertial
%   parameters EM learning
%  
%
% Author: Francesco Nori
% Genova, Dec 2014

% obj = updateNet(obj);    % called in the overload of setState
% obj = updateEvd(obj);    % called in the overload of setY


NB  = obj.IDmodel.modelParams.NB;

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

if strcmp(obj.inf_engine, 'jtree_inf_engine')
   %% jtree_inf_engine
   C = cliques_from_engine(obj.bnt.engine);
   I = cell2mat(obj.bnt.nodes.index);
   % J is the inverse index for the permutation I
   J(I) = 1:length(I);
   for i = 1:length(C)
      for j = 1 : length(C{i})
         if J(C{i}(j))>6*NB
            C{i}(j) = [];
         end
      end
      tmp       = marginal_nodes(obj.bnt.engine, C{i});
      obj.Sd_sm = set(obj.Sd_sm, tmp.Sigma, J(C{i}), J(C{i}));
   end
   ii = zeros(6*NB,1);
   for i = 1 : NB
      ii(     (i-1)*4+1:     4*i, 1) = (i-1)*6 + 1 : (i-1)*6 + 4;
      ii(4*NB+(i-1)*2+1:4*NB+2*i, 1) = (i-1)*6 + 5 : (i-1)*6 + 6;
   end
   obj.Sd = obj.Sd_sm(ii,ii);
   
   for i = 1:NB
      iSd(       (i-1)*4+1 :        4*i, 1) = [6 6 6 obj.IDmodel.jn(i)]';
      iSd(4*NB + (i-1)*2+1 : 4*NB + 2*i, 1) = [6 obj.IDmodel.jn(i)]';
   end
   
   Sd_sm = submatrix(iSd, iSd, obj.Sd);
   for i = 1 : NB
      S_ind((i-1)*6+1:(i-1)*6+4, 1) =        (i-1)*4+1:       i*4;
      S_ind((i-1)*6+5:(i-1)*6+6, 1) = 4*NB + (i-1)*2+1:4*NB + i*2;
   end
   obj.Sd = Sd_sm(S_ind, S_ind);
  
   
elseif strcmp(obj.inf_engine, 'gaussian_inf_engine')
   %% Gaussian inference engine
   I        = cell2mat(obj.bnt.nodes.index);
   ns       = obj.bnt.bnet.node_sizes;
   S        = marginal_nodes(obj.bnt.engine, I(1:6*NB));
   [~, p]   = sort(I(1:6*NB));
   ns       = ns(I(1:6*NB));
   ns       = ns(p);
   p_inv(p) = 1:length(p);
   
   obj.Sd_sm = submatrix(ns, ns, S.Sigma);
   obj.Sd    = obj.Sd_sm(p_inv,p_inv);

%    obj.Sd_sm = submatrix(obj.iSd, obj.iSd, obj.Sd);
%    for i = 1 : NB
%       S_ind(       (i-1)*4+1:       i*4, 1) = (i-1)*6+1:(i-1)*6+4;
%       S_ind(4*NB + (i-1)*2+1:4*NB + i*2, 1) = (i-1)*6+5:(i-1)*6+6;
%    end
%    obj.Sd = obj.Sd_sm(S_ind, S_ind);
%    
%    obj.Sd_sm = submatrix(obj.iSd, obj.iSd, obj.Sd);
%    for i = 1 : NB
%       S_ind((i-1)*6+1:(i-1)*6+4, 1) =        (i-1)*4+1:       i*4;
%       S_ind((i-1)*6+5:(i-1)*6+6, 1) = 4*NB + (i-1)*2+1:4*NB + i*2;
%    end
%    obj.Sd = obj.Sd_sm(S_ind, S_ind);

end

%% General computations
for i = 1:NB
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(1));
   a(i,:)   = tmp.mu';
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(2));
   fB(i,:)  = tmp.mu';
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(3));
   f(i,:)   = tmp.mu';
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(4));
   tau(i,1) = tmp.mu';
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(5));
   fx(i,:)  = tmp.mu';
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(6));
   d2q(i,1) = tmp.mu';
end



for i = 1 : NB
   obj.d((1:26)+(i-1)*26, 1) = [a( i,1:6), fB(i,1:6), f( i,1:6), tau(i,1), fx(i,1:6), d2q(i,1)]';
end

end % solveID