function obj = solveID(obj)
%solveID Differential Inverse Dynamics with Newton-Euler Algorithm (DNEA)
%   This function solves the differential inverse dynamics problem with the
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
%   link-i and d2q_i is acceleration of joint-i. The output x is structured
%   as follows:
%
%   x   = [q, dq]
%
%   where q are joint positions and dq joint accelerations. The input to 
%   the algorithm is in obj.IDmeas.y organized as follows:
%
%   obj.IDmeas.y = [y_1, y_2, ... , y_obj.IDsens.m]
%
%   The relationship between d and y is given by Y d = y where the
%   matrix Y is assumed constant with respect to q and dq. Moreover, the
%   variables d should satisfy the Newton-Euler equations represented as
%   D(x) d + b(x) = 0. The estimation of d and x is obtained by
%   a minimum variance estimation given the following definitions:
%
%   D(x_bar) d + dDb(x_bar, d_bar) x = 0
%
%          Y d - y                   = 0
%
%   where x_bar and d_bar is the value around which we linearize the
%   non-linear equations. 
%
% Author: Francesco Nori
% Genova, Dec 2014

NB = obj.IDmodel.modelParams.NB;
b  = sparse(obj.ibs, ones(length(obj.ibs),1), obj.bs, 19*NB, 1);
D  = sparse(obj.iDs, obj.jDs, obj.Ds, 19*NB, 26*NB); 

Dx = D(1:19*NB, 1:19*NB);
Dy = D(1:19*NB, 19*NB+1:26*NB);

% We write the estimation problem as:
%
%     [Dx Dy dDb] * [dx; dy; x] = v
%
% z - [Yx Yy   0] * [dx; dy; x] = w
%
%     x     ~ N(             mux, Sx);
%    dy     ~ N(             muy, Sy);
%     v     ~ N(               0, Sv);
% z|dx,dy,x ~ N([Yx Yy]*[dx; dy], Sw);
%
% With easy substitutions:
%
%     dx + Dx^(-1)*[Dy*dy + dDb*x] = Dx^(-1)v
%
%             z - [Yx Yy]*[dx; dy] = w
%
%         x ~ N(                     mux, Sx);
%        dy ~ N(                     muy, Sy);
%   dx|dy,x ~ N( Dx^(-1) [Dy dy + dDb x], Dx^(-1) Sv Dx^(-1)');
% z|dx,dy,x ~ N(   [Yx Yy 0]*[dx; dy; x], Sw);
%
% which is totally equivalent to 'SIXTH  EXAMPLE'
% in gaussSumImplicit.m Some simplifications using
% the Shur complement inversion have been used to
% reduce the computational cost.
%
%  S      = [Sy, -Sy*Dy'; -Dy*Sy, Sv + Dy*Sy*Dy'];
%  S^(-1) = [Sy^(-1)+Dy'*Sv^(-1)*Dy, Sv^(-1)*Dy; Dy'*Sv^(-1), Sv^(-1)];
% (S^(-1)+[Yy Yx]'*Sw^(-1)*[Yy Yx])^(-1)
%
% In practice computations are as before with the 
% following redefinitions:
%
%      Sy -> [Sy, 0; 0 Sx]
% [Yx Yy] -> [Yx Yy 0]
%      Dy -> [Dy dDb]

% Sv_inv = eye(19*NB)./sModel;
Sv_inv = obj.IDmodel.modelParams.Sv_inv;
% Sw_inv = eye(7*NB) ./sUknown;
Sw_inv = obj.IDmodel.modelParams.Sw_inv;
% Sy_inv = eye(my)   ./sMeas;
Sy_inv = obj.IDsens.sensorsParams.Sy_inv;

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
obj.d     = d;
% obj.Sd = full(inv(S1'*Ss*S1));
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