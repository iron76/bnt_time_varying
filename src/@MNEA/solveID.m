function obj = solveID(obj)
%solveID Inverse Dynamics with redundant measures Newton-Euler (MNEA)
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

% INIT: NON-SPARSE VERSION

% b  = obj.b(1:4*NB,1);
% iDx = zeros(4*NB,1);  iDy = zeros(2*NB,1);
% for i = 1:obj.IDmodel.modelParams.NB
%    iDx((i-1)*4+1 : i*4,1) = ((i-1)*6+1 : (i-1)*6+4)';
%    iDy((i-1)*2+1 : i*2,1) = ((i-1)*6+5 : (i-1)*6+6)';
% end
%
% Dx = obj.D(1:4*NB,iDx);
% Dy = obj.D(1:4*NB,iDy);

% INIT: SPARSE VERSION

D = sparse(obj.iDs, obj.jDs, obj.Ds, 19*NB, 26*NB);
b = sparse(obj.ibs, ones(size(obj.ibs)), obj.bs, 19*NB, 1);

Dx = D(1:19*NB, 1:19*NB);
Dy = D(1:19*NB, 19*NB+1:26*NB);

% Sv_inv = eye(19*NB)./sModel;
Sv_inv = obj.IDmodel.modelParams.Sv_inv.matrix;
% Sw_inv = eye(7*NB) ./sUknown;
Sw_inv = obj.IDmodel.modelParams.Sw_inv.matrix;
% Sy_inv = eye(my)   ./sMeas;
Sy_inv = obj.IDsens.sensorsParams.Sy_inv.matrix;

Sinv   = [Dx'*Sv_inv*Dx Dx'*Sv_inv*Dy; Dy'*Sv_inv*Dx, Sw_inv+ Dy'*Sv_inv*Dy];
Sw     = Sw_inv\sparse(1:7*NB, 1:7*NB, 1);
Dx_inv = Dx\sparse(1:19*NB, 1:19*NB, 1);
Y = obj.IDsens.sensorsParams.Ys;

Ss = Sinv+Y'*Sy_inv*Y;
% L = chol(S1'*Ss*S1, 'lower');    % S1'*W*S1 = L*L'
% PWinv1 = S1*inv_chol(L)*S1';
% Ls = S1*L;

% Sxy = [Sv + Dy*Sy*Dy', -Dy*Sy; -Sy*Dy', Sy];
% mxy = -Sxy*[ -Sv^(-1)*mx; -Dy'*Sv^(-1)*mx - Sy^(-1)*my]


% Sxy = S = [Dx^(-1)*inv(Sv_inv)*Dx^(-1)' + Dx^(-1)*Dy*inv(Sw_inv)*Dy'*Dx^(-1)', -Dx^(-1)*Dy*inv(Sw_inv); -inv(Sw_inv)*Dy'*Dx^(-1)', inv(Sw_inv)],1)
Sxy = [Dx_inv + Dx_inv*Dy*Sw*Dy'*Sv_inv, -Dx_inv*Dy*Sw; -Sw*Dy'*Sv_inv, Sw];
mx  = -b;
my  = zeros(7*NB,1);
mxy = -Sxy*[-mx; -Dy'*Sv_inv*mx - Sw_inv*my];
% d   = mxy + Ss\Y'*Sy_inv*(obj.IDmeas.y-Y*mxy);
[~,~,S1] = chol(Ss, 'lower');
d   = mxy +S1*((S1'*Ss*S1)\(S1'*(Y'*Sy_inv*(obj.IDmeas.y-Y*mxy))));

dx    =      d(1:NB*19      , 1);
dy    =      d(1+NB*19 : end, 1);

dxc  = mat2cell( dx, 19*ones(1, NB), 1);
dyc  = mat2cell( dy,  7*ones(1, NB), 1);

for i = 1 : NB
   dc{i}   = [dxc{i,1}; dyc{i,1}];
   
   obj.a  (1:6,i) = dc{i}( 1: 6, 1);
   obj.fB (1:6,i) = dc{i}( 7: 12, 1);
   obj.f  (1:6,i) = dc{i}(13: 18, 1);
   obj.tau(1:1,i)  = dc{i}(19, 1);
   obj.fx (1:6,i) = dc{i}( 20: 25, 1);
   obj.d2q(1:1,i)  = dc{i}(26, 1);
   
   d((1:26)+(i-1)*26, 1) = [obj.a(1:6,i); obj.fB(1:6,i); obj.f(1:6,i); obj.tau(1,i); obj.fx(1:6,i); obj.d2q(1,i)];
end

obj.d  = d;
obj.Sd = full(inv(Ss));

for i = 1:NB
   iSd(       (i-1)*4+1 :        4*i, 1) = [6 6 6 obj.IDmodel.jn(i)]';
   iSd(4*NB + (i-1)*2+1 : 4*NB + 2*i, 1) = [6 obj.IDmodel.jn(i)]';
end

obj.Sd_sm = submatrix(iSd, iSd, obj.Sd);
for i = 1 : NB
   S_ind((i-1)*6+1:(i-1)*6+4, 1) =        (i-1)*4+1:       i*4;
   S_ind((i-1)*6+5:(i-1)*6+6, 1) = 4*NB + (i-1)*2+1:4*NB + i*2;
end
obj.Sd = obj.Sd_sm(S_ind, S_ind);

end % solveID