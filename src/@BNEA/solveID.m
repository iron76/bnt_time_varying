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

for i = NB:-1:1
   tmp      = marginal_nodes(obj.bnt.engine, obj.bnt.nodes.index{i}(1));
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