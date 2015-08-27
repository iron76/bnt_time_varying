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

% Dx = D(1:19*NB, 1:19*NB);
% Dy = D(1:19*NB, 19*NB+1:26*NB);

% Sv_inv = eye(19*NB)./sModel;
Sv_inv = obj.IDmodel.modelParams.Sv_inv.matrix;
% Sw_inv = eye(7*NB) ./sUknown;
Sw_inv = obj.IDmodel.modelParams.Sw_inv.matrix;
% Sw     = obj.IDmodel.modelParams.Sw.matrix;
% Sy_inv = eye(my)   ./sMeas;
Sy_inv = obj.IDsens.sensorsParams.Sy_inv.matrix;

Y = obj.IDsens.sensorsParams.Ys;

y      = obj.IDmeas.y;
S_Dinv = Sv_inv;
S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);
S_Yinv = Sy_inv;
bY     = zeros(size(y));
bD     = b;
muD    = zeros(length(S_dinv), 1);
d      = (D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y)\(Y'*S_Yinv*(y-bY) - D'*S_Dinv*bD + S_dinv * muD);

% shuffle from [dx dy] to d
obj.d  = d(obj.id,1);
obj.Sd = full(inv(D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y));

% shuffle from Sd from [dx dy] to d
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