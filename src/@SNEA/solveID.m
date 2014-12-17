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

%%
%
ib  = zeros(obj.IDmodel.sparseParams.nb, 1);
b1s = zeros(obj.IDmodel.sparseParams.nb, 1);

iD  = zeros(sum(obj.IDmodel.sparseParams.nD), 1);
jD  = zeros(sum(obj.IDmodel.sparseParams.nD), 1);
D1s = zeros(sum(obj.IDmodel.sparseParams.nD), 1);

pD = 1;

for i = 1:obj.IDmodel.modelParams.NB  
  if obj.IDmodel.modelParams.parent(i) == 0
    % a{i} = obj.Xup{i}*(-a_grav) + obj.IDmodel.S{i}*qdd(i);
    % D1  = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, obj.jn(i)) zeros(6,6) obj.IDmodel.S{i}];
    % b1  = obj.Xup{i}*(-a_grav);
    b1s((i-1)*12+1: (i-1)*12+6, 1) = obj.Xup{i}*(-obj.IDmodel.g);
    ib((i-1)*12+1: (i-1)*12+6, 1)  = (i-1)*19+1: (i-1)*19+6;

    D1s(pD : pD+5, 1) = -1*ones(6,1);
    iD (pD : pD+5, 1) = obj.IDmodel.sparseParams.iD11(1:6:36,i,i)+[0 1 2 3 4 5]';
    jD (pD : pD+5, 1) = obj.IDmodel.sparseParams.jD11(1:6:36,i,i);
    pD = pD + 6;
    
    D1s(pD: pD+obj.IDmodel.jn(i)*6-1, 1) = obj.IDmodel.S{i};
    iD (pD: pD+obj.IDmodel.jn(i)*6-1, 1) = obj.IDmodel.sparseParams.iD16(:,i,i);
    jD (pD: pD+obj.IDmodel.jn(i)*6-1, 1) = obj.IDmodel.sparseParams.jD16(:,i,i);  
    pD = pD + obj.IDmodel.jn(i)*6;
  else
    % a{i} = ... + obj.IDmodel.S{i}*qdd(i) + crm(obj.v(:,i))*vJ;
    % vJ = obj.IDmodel.S{i}*qd(i);
    % D1 = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, obj.jn(i)) zeros(6,6) obj.IDmodel.S{i}];
    % b1 = crm(obj.v(:,i))*vJ;
    
    b1s((i-1)*12+1: (i-1)*12+6, 1)    = crm(obj.v(:,i))*obj.vJ(:,i);
    ib((i-1)*12+1: (i-1)*12+6, 1)     = (i-1)*19+1: (i-1)*19+6;
    
    D1s(pD : pD+5, 1) = -1*ones(6,1);
    iD (pD : pD+5, 1) = obj.IDmodel.sparseParams.iD11(1:6:36,i,i)+[0 1 2 3 4 5]';
    jD (pD : pD+5, 1) = obj.IDmodel.sparseParams.jD11(1:6:36,i,i);
    pD = pD + 6;
    
    D1s(pD: pD+obj.IDmodel.jn(i)*6-1, 1) = obj.IDmodel.S{i};
    iD (pD: pD+obj.IDmodel.jn(i)*6-1, 1) = obj.IDmodel.sparseParams.iD16(:,i,i);
    jD (pD: pD+obj.IDmodel.jn(i)*6-1, 1) = obj.IDmodel.sparseParams.jD16(:,i,i);  
    pD = pD + obj.IDmodel.jn(i)*6;
    % a{i} = obj.Xup{i}*a{obj.IDmodel.modelParams.parent(i)} + ...
    % Dc{i, obj.IDmodel.modelParams.parent(i)} = [ obj.Xup{i} zeros(6,6) zeros(6,6) zeros(6, obj.jn(i)) zeros(6,6) zeros(6, obj.jn(i))
    %     zeros(12+obj.jn(i), 24+2*obj.jn(i))];
    
    j = obj.IDmodel.modelParams.parent(i);
    D1s(pD: pD+35, 1) = obj.Xup{i}(:);
    iD (pD: pD+35, 1) = obj.IDmodel.sparseParams.iD11(:,j,i);
    jD (pD: pD+35, 1) = obj.IDmodel.sparseParams.jD11(:,j,i);  
    pD = pD + 36;    

  end
  % fB{i} = obj.IDmodel.modelParams.I{i}*a{i} + crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i);
  % D2 = [obj.IDmodel.modelParams.I{i} -eye(6) zeros(6,6) zeros(6, obj.jn(i)) zeros(6,6) zeros(6, obj.jn(i))];
  % b2 = crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i);
  
  b1s((i-1)*12+7: i*12, 1) = crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i);
  ib((i-1)*12+7: i*12, 1)  = (i-1)*19+7: (i-1)*19+12;

  D1s(pD : pD+35, 1) = obj.IDmodel.modelParams.I{i}(:);
  iD (pD : pD+35, 1) = obj.IDmodel.sparseParams.iD21(:,i,i);
  jD (pD : pD+35, 1) = obj.IDmodel.sparseParams.jD21(:,i,i);
  pD = pD + 36;
  
  D1s(pD : pD+5, 1) = -1*ones(6,1);
  iD (pD : pD+5, 1) = obj.IDmodel.sparseParams.iD22(1:6:36,i,i)+[0 1 2 3 4 5]';
  jD (pD : pD+5, 1) = obj.IDmodel.sparseParams.jD22(1:6:36,i,i);
  pD = pD + 6;
  
  % f{i} = fB{i} - obj.Xa{i}' \ f_ext{i};
  % f{obj.IDmodel.modelParams.parent(j)} = f{obj.IDmodel.modelParams.parent(j)} + obj.Xup{j}'*f{j};
  % D3 = [zeros(6,6) eye(6) -eye(6) zeros(6, obj.jn(i)) -inv(obj.Xa{i}') zeros(6, obj.jn(i))];
  % b3 = zeros(6,1);
  
  A  = -inv(obj.Xa{i}');

  D1s(pD : pD+5, 1) = 1*ones(6,1);
  iD (pD : pD+5, 1) = obj.IDmodel.sparseParams.iD32(1:6:36,i,i)+[0 1 2 3 4 5]';
  jD (pD : pD+5, 1) = obj.IDmodel.sparseParams.jD32(1:6:36,i,i);
  pD = pD + 6;

  D1s(pD : pD+5, 1) = -1*ones(6,1);
  iD (pD : pD+5, 1) = obj.IDmodel.sparseParams.iD33(1:6:36,i,i)+[0 1 2 3 4 5]';
  jD (pD : pD+5, 1) = obj.IDmodel.sparseParams.jD33(1:6:36,i,i);
  pD = pD + 6;

  D1s(pD : pD+35, 1) = A(:);
  iD (pD : pD+35, 1) = obj.IDmodel.sparseParams.iD35(:,i,i);
  jD (pD : pD+35, 1) = obj.IDmodel.sparseParams.jD35(:,i,i);
  pD = pD + 36;
    
  % tau(i,1) = obj.IDmodel.S{i}' * f{i};
  % D4 = [zeros(obj.jn(i), 6) zeros(obj.jn(i), 6) obj.IDmodel.S{i}' -eye(obj.jn(i)) zeros(obj.jn(i), 6) zeros(obj.jn(i), obj.jn(i))];
  % b4 =  zeros(obj.jn(i), 1);
  
  D1s(pD : pD+5, 1) = obj.IDmodel.S{i}';
  iD (pD : pD+5, 1) = obj.IDmodel.sparseParams.iD43(:,i,i);
  jD (pD : pD+5, 1) = obj.IDmodel.sparseParams.jD43(:,i,i);
  pD = pD + 6*obj.IDmodel.jn(i);
  
  D1s(pD : pD+obj.IDmodel.jn(i)-1, 1) = -ones(1,obj.IDmodel.jn(i));
  iD (pD : pD+obj.IDmodel.jn(i)-1, 1) = obj.IDmodel.sparseParams.iD44(:,i,i);
  jD (pD : pD+obj.IDmodel.jn(i)-1, 1) = obj.IDmodel.sparseParams.jD44(:,i,i);
  pD = pD + obj.IDmodel.jn(i);

  for j = obj.IDmodel.sparseParams.ind_j{i}
    % f{obj.IDmodel.modelParams.parent(j)} = f{obj.IDmodel.modelParams.parent(j)} + obj.Xup{j}'*f{j};
    % Dc{i,j} = [ zeros(12, 24+2*obj.jn(i))
    %     zeros(6,6) zeros(6,6) obj.Xup{j}' zeros(6, obj.jn(i)) zeros(6,6) zeros(6, obj.jn(i))
    %     zeros(obj.jn(i), 24+2*obj.jn(i))];
    
    A       = obj.Xup{j}';
    
    D1s(pD : pD+35, 1) = A(:);
    iD (pD : pD+35, 1) = obj.IDmodel.sparseParams.iD33(:,j,i);
    jD (pD : pD+35, 1) = obj.IDmodel.sparseParams.jD33(:,j,i);
    pD = pD + 36;
  end
end

NB = obj.IDmodel.modelParams.NB;
b  = sparse(ib, ones(length(ib),1), b1s, 19*NB, 1);
Ds = sparse(iD, jD, D1s, 19*NB, 26*NB); 

Dx = Ds(1:19*NB, 1:19*NB);
Dy = Ds(1:19*NB, 19*NB+1:26*NB);

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
d   = mxy + Ss\Y'*Sy_inv*(obj.IDmeas.y-Y*mxy);
% d   = mxy +S1*((S1'*Ss*S1)\(S1'*(Y'*Sy_inv*(y-Y*mxy))));

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
obj.d = d;


end % solveID