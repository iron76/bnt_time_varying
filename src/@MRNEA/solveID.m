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


NB = obj.IDmodel.modelParams.NB;
b  = sparse(obj.ibs, ones(length(obj.ibs),1), obj.bs, 19*NB, 1);
D  = sparse(obj.iDs, obj.jDs, obj.Ds, 19*NB, 26*NB); 

Dx = D(1:19*NB, 1:19*NB);
Dy = D(1:19*NB, 19*NB+1:26*NB);

% Dx_ind = zeros(4*NB,1);
% Dy_ind = zeros(2*NB,1);
% for i = 1 : NB
%    Dx_ind((i-1)*4+1:i*4, 1) = (i-1)*6+1  : (i-1)*6+4;
%    Dy_ind((i-1)*2+1:i*2, 1) = (i-1)*6+5  : (i-1)*6+6;
% end
% 
% size(obj.D(1:4*NB,1:6*NB))
% norm([obj.D(1:4*NB,Dx_ind) obj.D(1:4*NB,Dy_ind)] - [ Dx Dy ])

dy = obj.IDmeas.y;

iy = 1;
for j = 1 : obj.IDsens.sensorsParams.ny
   for i = 1 : obj.IDmodel.n
      if strcmp(obj.IDsens.sensorsParams.labels{j,1}, ['fx' num2str(i)])
         dy(7*(i-1)+1 : 7*(i-1)+6,1)  = obj.IDmeas.y(iy:iy+5,1);
      elseif strcmp(obj.IDsens.sensorsParams.labels{j,1}, ['d2q' num2str(i)])
         dy(7*i,1) = obj.IDmeas.y(iy, 1);
      end
   end
   iy = iy + obj.IDsens.sensorsParams.sizes{j,1};
end

dx = -Dx\(Dy*dy+b);
for i = 1 : NB
   obj.d((i-1)*26+1  : (i-1)*26+19, 1) = dx((i-1)*19+1:i*19,1);
   obj.d((i-1)*26+20 : i*26, 1)     = dy((i-1)*7+1:i*7,1);
end % solveID