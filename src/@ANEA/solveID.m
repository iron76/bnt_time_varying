function obj = solveID(obj)
%solveID accelration with sparse Newton-Euler Algorithm (SNEA)
%   This function solves the accelerations estimation with the sparse
%   Newton-Euler algorithm, as described in the paper "BERDY: Bayesian 
%   Estimation for Robot Dynamics. A Probabilistic Estimation of Whole-Body
%   Dynamics with Redundant Measurements." The output 'd' is structured as
%   follows:
%
%   d   = [d_1, d_2, ..., d_obj.IDstate.n]
%
%   where:
%
%   d_i = [a_i, d2q_i]
%
%   and a_i is the link-i spatial accelration, d2q_i is acceleration of 
%   joint-i. The input to the algorithm is in obj.IDmeas.y organized as 
%   follows:
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

D = obj.D.matrix;
b = obj.b.matrix;
Sv_inv = obj.IDmodel.modelParams.Sv_inv.matrix;
Sw_inv = obj.IDmodel.modelParams.Sw_inv.matrix;
Sy_inv = obj.IDsens.sensorsParams.Sy_inv.matrix;
Y      = obj.IDsens.sensorsParams.Ys;
Y      = Y(1:obj.IDmeas.m, 1:7*NB);
Y      = Y(:, obj.id);

S_Dinv = Sv_inv;
S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);
S_Yinv = Sy_inv;
bY     = zeros(obj.IDmeas.m,1);
bD     = b;
S_dinv = S_dinv(obj.id,obj.id);

obj.Sd = inv(D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y);
obj.d  = (D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y)\(Y'*S_Yinv*(obj.IDmeas.y-bY) - D'*S_Dinv*bD);

end % solveID