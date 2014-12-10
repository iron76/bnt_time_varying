function  [d, S1] = sID( model, q, qd , Yx, Yy, y, S1, sparseModel)

% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

% n  = model.NB;
% Size of the D matrix [dm, dn]
% dm = 18*model.NB + n
% dn = 24 * model.NB + 2*n

global dSv_inv idSv_inv jdSv_inv
global dSw_inv idSw_inv jdSw_inv
global dSy_inv idSy_inv jdSy_inv

NB  = model.NB;

a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  [~, jn{i}] = size(S{i});
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    % a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    % a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
  end
end

if nargin >= 5
  [~, Xa] = apply_external_forces( model.parent, Xup, cell(NB), cell(NB) );
end


ib  = zeros(sparseModel.nb, 1);
b1s = zeros(sparseModel.nb, 1);

iD  = zeros(sum(sparseModel.nD), 1);
jD  = zeros(sum(sparseModel.nD), 1);
D1s = zeros(sum(sparseModel.nD), 1);

pD = 1;

for i = 1:model.NB  
  if model.parent(i) == 0
    % a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
    % D1  = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
    % b1  = Xup{i}*(-a_grav);
    b1s((i-1)*12+1: (i-1)*12+6, 1) = Xup{i}*(-a_grav);
    ib((i-1)*12+1: (i-1)*12+6, 1)  = (i-1)*19+1: (i-1)*19+6;

    D1s(pD : pD+5, 1) = -1*ones(6,1);
    iD (pD : pD+5, 1) = sparseModel.iD11(1:6:36,i,i)+[0 1 2 3 4 5]';
    jD (pD : pD+5, 1) = sparseModel.jD11(1:6:36,i,i);
    pD = pD + 6;
    
    D1s(pD: pD+jn{i}*6-1, 1) = S{i}(:);
    iD (pD: pD+jn{i}*6-1, 1) = sparseModel.iD16(:,i,i);
    jD (pD: pD+jn{i}*6-1, 1) = sparseModel.jD16(:,i,i);  
    pD = pD + jn{i}*6;
  else
    % a{i} = ... + S{i}*qdd(i) + crm(v{i})*vJ;
    vJ = S{i}*qd(i);
    % D1 = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
    % b1 = crm(v{i})*vJ;
    
    b1s((i-1)*12+1: (i-1)*12+6, 1)    = crm(v{i})*vJ;
    ib((i-1)*12+1: (i-1)*12+6, 1)     = (i-1)*19+1: (i-1)*19+6;
    
    D1s(pD : pD+5, 1) = -1*ones(6,1);
    iD (pD : pD+5, 1) = sparseModel.iD11(1:6:36,i,i)+[0 1 2 3 4 5]';
    jD (pD : pD+5, 1) = sparseModel.jD11(1:6:36,i,i);
    pD = pD + 6;
    
    D1s(pD: pD+jn{i}*6-1, 1) = S{i}(:);
    iD (pD: pD+jn{i}*6-1, 1) = sparseModel.iD16(:,i,i);
    jD (pD: pD+jn{i}*6-1, 1) = sparseModel.jD16(:,i,i);  
    pD = pD + jn{i}*6;
    % a{i} = Xup{i}*a{model.parent(i)} + ...
    % Dc{i, model.parent(i)} = [ Xup{i} zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})
    %     zeros(12+jn{i}, 24+2*jn{i})];
    
    j = model.parent(i);
    D1s(pD: pD+35, 1) = Xup{i}(:);
    iD (pD: pD+35, 1) = sparseModel.iD11(:,j,i);
    jD (pD: pD+35, 1) = sparseModel.jD11(:,j,i);  
    pD = pD + 36;    

  end
  % fB{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
  % D2 = [model.I{i} -eye(6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})];
  % b2 = crf(v{i})*model.I{i}*v{i};
  
  b1s((i-1)*12+7: i*12, 1) = crf(v{i})*model.I{i}*v{i};
  ib((i-1)*12+7: i*12, 1)  = (i-1)*19+7: (i-1)*19+12;

  D1s(pD : pD+35, 1) = model.I{i}(:);
  iD (pD : pD+35, 1) = sparseModel.iD21(:,i,i);
  jD (pD : pD+35, 1) = sparseModel.jD21(:,i,i);
  pD = pD + 36;
  
  D1s(pD : pD+5, 1) = -1*ones(6,1);
  iD (pD : pD+5, 1) = sparseModel.iD22(1:6:36,i,i)+[0 1 2 3 4 5]';
  jD (pD : pD+5, 1) = sparseModel.jD22(1:6:36,i,i);
  pD = pD + 6;
  
  % f{i} = fB{i} - Xa{i}' \ f_ext{i};
  % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};
  % D3 = [zeros(6,6) eye(6) -eye(6) zeros(6, jn{i}) -inv(Xa{i}') zeros(6, jn{i})];
  % b3 = zeros(6,1);
  
  A  = -inv(Xa{i}');

  D1s(pD : pD+5, 1) = 1*ones(6,1);
  iD (pD : pD+5, 1) = sparseModel.iD32(1:6:36,i,i)+[0 1 2 3 4 5]';
  jD (pD : pD+5, 1) = sparseModel.jD32(1:6:36,i,i);
  pD = pD + 6;

  D1s(pD : pD+5, 1) = -1*ones(6,1);
  iD (pD : pD+5, 1) = sparseModel.iD33(1:6:36,i,i)+[0 1 2 3 4 5]';
  jD (pD : pD+5, 1) = sparseModel.jD33(1:6:36,i,i);
  pD = pD + 6;

  D1s(pD : pD+35, 1) = A(:);
  iD (pD : pD+35, 1) = sparseModel.iD35(:,i,i);
  jD (pD : pD+35, 1) = sparseModel.jD35(:,i,i);
  pD = pD + 36;
    
  % tau(i,1) = S{i}' * f{i};
  % D4 = [zeros(jn{i}, 6) zeros(jn{i}, 6) S{i}' -eye(jn{i}) zeros(jn{i}, 6) zeros(jn{i}, jn{i})];
  % b4 =  zeros(jn{i}, 1);
  
  D1s(pD : pD+5, 1) = S{i}';
  iD (pD : pD+5, 1) = sparseModel.iD43(:,i,i);
  jD (pD : pD+5, 1) = sparseModel.jD43(:,i,i);
  pD = pD + 6*jn{i};
  
  D1s(pD : pD+jn{i}-1, 1) = -ones(1,jn{i});
  iD (pD : pD+jn{i}-1, 1) = sparseModel.iD44(:,i,i);
  jD (pD : pD+jn{i}-1, 1) = sparseModel.jD44(:,i,i);
  pD = pD + jn{i};

  for j = sparseModel.ind_j{i}
    % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};
    % Dc{i,j} = [ zeros(12, 24+2*jn{i})
    %     zeros(6,6) zeros(6,6) Xup{j}' zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})
    %     zeros(jn{i}, 24+2*jn{i})];
    
    A       = Xup{j}';
    
    D1s(pD : pD+35, 1) = A(:);
    iD (pD : pD+35, 1) = sparseModel.iD33(:,j,i);
    jD (pD : pD+35, 1) = sparseModel.jD33(:,j,i);
    pD = pD + 36;
  end
end

b  = sparse(ib, ones(length(ib),1), b1s, 19*NB, 1);
Ds = sparse(iD, jD, D1s, 19*NB, 26*NB); 

[my, ~]   = size(y);

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
Sv_inv = sparse(idSv_inv, jdSv_inv, dSv_inv);
% Sy_inv = eye(7*NB) ./sUknown;
Sy_inv = sparse(idSw_inv, jdSw_inv, dSw_inv);
% Sw_inv = eye(my)   ./sMeas;
Sw_inv = sparse(idSy_inv, jdSy_inv, dSy_inv);
Sinv   = [Dx'*Sv_inv*Dx Dx'*Sv_inv*Dy; Dy'*Sv_inv*Dx, Sy_inv+ Dy'*Sv_inv*Dy];

Dx_inv = Dx\sparse(1:19*NB, 1:19*NB, 1);
Sy     = Sy_inv\sparse(1:7*NB, 1:7*NB, 1);

Ss = Sinv+[Yx Yy]'*Sw_inv*[Yx Yy];
% L = chol(S1'*Ss*S1, 'lower');    % S1'*W*S1 = L*L'
% PWinv1 = S1*inv_chol(L)*S1';
% Ls = S1*L;


% Sxy = [Sv + Dy*Sy*Dy', -Dy*Sy; -Sy*Dy', Sy];
% mxy = -Sxy*[ -Sv^(-1)*mx; -Dy'*Sv^(-1)*mx - Sy^(-1)*my]


% Sxy = S = [Dx^(-1)*inv(Sv_inv)*Dx^(-1)' + Dx^(-1)*Dy*inv(Sy_inv)*Dy'*Dx^(-1)', -Dx^(-1)*Dy*inv(Sy_inv); -inv(Sy_inv)*Dy'*Dx^(-1)', inv(Sy_inv)],1)
Sxy = [Dx_inv + Dx_inv*Dy*Sy*Dy'*Sv_inv, -Dx_inv*Dy*Sy; -Sy*Dy'*Sv_inv, Sy];
mx  = -b;
my  = zeros(7*NB,1);
mxy = -Sxy*[-mx; -Dy'*Sv_inv*mx - Sy_inv*my];
% d   = mxy + Ss\[Yx Yy]'*Sw_inv*(y-[Yx Yy]*mxy);
d   = mxy +S1*((S1'*Ss*S1)\(S1'*([Yx Yy]'*Sw_inv*(y-[Yx Yy]*mxy))));

dx    =      d(1:NB*19      , 1);
dy    =      d(1+NB*19 : end, 1);

dxc  = mat2cell( dx, 19*ones(1, NB), 1);
dyc  = mat2cell( dy,  7*ones(1, NB), 1);

a   = zeros(NB,6);
fB  = zeros(NB,6);
f   = zeros(NB,6);
tau = zeros(NB,1);
d2q = zeros(NB,1);
fx  = zeros(NB,6);

for i = 1 : NB
  dc{i}   = [dxc{i,1}; dyc{i,1}];
  
  a( i,1:6) = dc{i}( 1: 6, 1);
  fB(i,1:6) = dc{i}( 7: 12, 1);
  f( i,1:6) = dc{i}(13: 18, 1);
  tau(i,1)  = dc{i}(19, 1);
  fx(i,1:6) = dc{i}( 20: 25, 1);
  d2q(i,1)  = dc{i}(26, 1);
  
  d((1:26)+(i-1)*26, 1) = [a( i,1:6), fB(i,1:6), f( i,1:6), tau(i,1), fx(i,1:6), d2q(i,1)]';
end
