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
b  = obj.b.matrix;

Dx = obj.D(1:NB, 1:2:2*NB);
Dy = obj.D(1:NB, 2:2:2*NB);

% We write the estimation problem as:
%
%     [Dx Dy]*[dx; dy] = v
%
% z - [Yx Yy]*[dx; dy] = w
%
%    dy   ~ N(             muy, Sy);
%     v   ~ N(               0, Sv);
% z|dx,dy ~ N([Yx Yy]*[dx; dy], Sw);
%
% With easy substitutions:
%
%     dx + Dx^(-1)*Dy*dy] = Dx^(-1)v
%
%    z - [Yx Yy]*[dx; dy] = w
%
%   dy    ~ N(             muy, Sy);
%dx|dy    ~ N(               0, Dx^(-1)*Sv*Dx^(-1)');
% z|dx,dy ~ N([Yx Yy]*[dx; dy], Sw);
%
% which is totally equivalent to 'SIXTH  EXAMPLE'
% in gaussSumImplicit.m Some simplifications using`
% the Shur complement inversion have been used to
% reduce the computational cost.
%
%  S      = [Sy, -Sy*Dy'; -Dy*Sy, Sv + Dy*Sy*Dy'];
%  S^(-1) = [Sy^(-1)+Dy'*Sv^(-1)*Dy, Sv^(-1)*Dy; Dy'*Sv^(-1), Sv^(-1)];
% (S^(-1)+[Yy Yx]'*Sw^(-1)*[Yy Yx])^(-1)

% Sv_inv = eye(19*NB)./sModel;
Sv_inv = obj.IDmodel.modelParams.Sv_inv.matrix;
% Sw_inv = eye(7*NB) ./sUknown;
Sw_inv = obj.IDmodel.modelParams.Sw_inv.matrix;
% Sw     = obj.IDmodel.modelParams.Sw.matrix;
% Sy_inv = eye(my)   ./sMeas;
if isa(obj.IDsens.sensorsParams.Sy_inv, 'submatrixSparse')
   Sy_inv = obj.IDsens.sensorsParams.Sy_inv.matrix;
else
   Sy_inv = obj.IDsens.sensorsParams.Sy_inv;
end

Sinv   = [Dx'*Sv_inv*Dx Dx'*Sv_inv*Dy; Dy'*Sv_inv*Dx, Sw_inv+ Dy'*Sv_inv*Dy];
% Dx_inv = Dx\sparse(1:19*NB, 1:19*NB, 1);
Y = obj.IDsens.sensorsParams.Ys;


Ss = Sinv+Y'*Sy_inv*Y;
% L = chol(S1'*Ss*S1, 'lower');    % S1'*W*S1 = L*L'
% PWinv1 = S1*inv_chol(L)*S1';
% Ls = S1*L;

% Sxy = [Sv + Dy*Sy*Dy', -Dy*Sy; -Sy*Dy', Sy];
% mxy = -Sxy*[ -Sv^(-1)*mx; -Dy'*Sv^(-1)*mx - Sy^(-1)*my]


% Sxy = S = [Dx^(-1)*inv(Sv_inv)*Dx^(-1)' + Dx^(-1)*Dy*inv(Sw_inv)*Dy'*Dx^(-1)', -Dx^(-1)*Dy*inv(Sw_inv); -inv(Sw_inv)*Dy'*Dx^(-1)', inv(Sw_inv)],1)
% Sxy = [Dx_inv + Dx_inv*Dy*Sw*Dy'*Sv_inv, -Dx_inv*Dy*Sw; -Sw*Dy'*Sv_inv, Sw];
mx  = -b;
my  = zeros(1*NB,1);
% mxy = -Sxy*[-mx; -Dy'*Sv_inv*mx - Sw_inv*my];
mxy = [Dx\(mx-Dy*my); my];
d   = mxy + Ss\Y'*Sy_inv*(obj.IDmeas.y-Y*mxy);
% [~,~,S1] = chol(Ss, 'lower');
% d   = mxy +S1*((S1'*Ss*S1)\(S1'*(Y'*Sy_inv*(obj.IDmeas.y-Y*mxy))));

% shuffle from [dx dy] to d
obj.d = d(obj.id,1);

% obj.Sd = full(inv(S1'*Ss*S1));
obj.Sd = full(inv(Ss));


for i = 1:NB
   iSd(     i, 1) = 6;
   iSd(NB + i, 1) = obj.IDmodel.jn(i);
end
obj.Sd_sm = submatrix(iSd, iSd, obj.Sd);
for i = 1 : NB
   S_ind((i-1)*2+1:(i-1)*2+1, 1) = i;
   S_ind((i-1)*2+2:(i-1)*2+2, 1) = NB + i;
end
obj.Sd = obj.Sd_sm(S_ind, S_ind);

end % solveID