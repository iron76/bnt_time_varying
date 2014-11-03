function  [tau, a, fB, f, Sf] = mID( model, q, qd , Yc , y)

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

global sModel sMeas sUknown

NB  = model.NB;

Dc = cell(model.NB, model.NB);
bc = cell(model.NB, model.NB);

a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  [~, jn{i}] = size(S{i});
  for j = 1:model.NB
    Dc{i,j} = zeros(18+jn{i}, 24+2*jn{i});
  end
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


for i = 1:model.NB
  if model.parent(i) == 0
    % a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
    D1  = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
    b1  = Xup{i}*(-a_grav);
  else
    % a{i} = ... + S{i}*qdd(i) + crm(v{i})*vJ;
    D1 = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
    vJ = S{i}*qd(i);
    b1 = crm(v{i})*vJ;
    
    % a{i} = Xup{i}*a{model.parent(i)} + ...
    Dc{i, model.parent(i)} = [ Xup{i} zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})
        zeros(12+jn{i}, 24+2*jn{i})];
  end
  % fB{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
  D2 = [model.I{i} -eye(6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})];
  b2 = crf(v{i})*model.I{i}*v{i};
  
  Dc{i,i} = [D1; D2];
  bc{i}   = [b1; b2];
end

for i = model.NB:-1:1
  % f{i} = fB{i} - Xa{i}' \ f_ext{i};
  % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};
  D3 = [zeros(6,6) eye(6) -eye(6) zeros(6, jn{i}) -inv(Xa{i}') zeros(6, jn{i})];
  b3 = zeros(6,1);
  ind_j  = find(model.parent == i);
  for j = ind_j
    % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};  
    Dc{i,j} = [ zeros(12, 24+2*jn{i})
        zeros(6,6) zeros(6,6) Xup{j}' zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})
        zeros(jn{i}, 24+2*jn{i})];
  end

  % tau(i,1) = S{i}' * f{i};
  D4 = [zeros(jn{i}, 6) zeros(jn{i}, 6) S{i}' -eye(jn{i}) zeros(jn{i}, 6) zeros(jn{i}, jn{i})];
  b4 =  zeros(jn{i}, 1);
  
  Dc{i,i} = [Dc{i,i}; D3; D4];
  bc{i}   = [bc{i};   b3; b4];
end

yc = y.values;
D  = cell2mat(Dc);
Y  = cell2mat(Yc);
b  = cell2mat(bc);
y  = cell2mat(yc);
% d  = pinv([D; Y])*[-b; y];
A  = [D; Y];
d  = (A'*A)\A'*[-b; y];

dc  = mat2cell(d, ones(NB,1).*26 , 1);
a   = zeros(NB,6);
fB  = zeros(NB,6);
f   = zeros(NB,6);
tau = zeros(NB,1);
d2q = zeros(NB,1);
fx  = zeros(NB,6);

for i = 1 : NB
  a( i,1:6) = dc{i,1}( 1: 6);
  fB(i,1:6) = dc{i,1}( 7: 12);
  f( i,1:6) = dc{i,1}(13: 18);  
  tau(i,1)  = dc{i,1}(19);
  fx(i,1:6) = dc{i,1}( 20: 25);
  d2q(i,1)  = dc{i,1}(26);
end

[ny, ~]   = size(Yc);
[my, ~]   = size(y);
Dxc = cell(NB,NB);
Dyc = cell(NB,NB);
Yxc = cell(ny,NB);
Yyc = cell(ny,NB);
for i = 1 : NB
  for j = 1 : NB
    Dxc{i,j} = Dc{i,j}(:, 1:19);
    Dyc{i,j} = Dc{i,j}(:, 20:end);
  end
end
for i = 1 : ny
  for j = 1 : NB
    Yxc{i,j}  = Yc{i,j}(:, 1:19);
    Yyc{i,j}  = Yc{i,j}(:, 20:end);
  end
end

Dx = cell2mat(Dxc);
Dy = cell2mat(Dyc);
Yx = cell2mat(Yxc);
Yy = cell2mat(Yyc);

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
% in gaussSumImplicit.m Some simplifications using
% the Shur complement inversion have been used to
% reduce the computational cost.
%
% S = [Sy, -Sy*Dy'; -Dy*Sy, Sv + Dy*Sy*Dy'];
%
% and the solution is:
%
% (S^(-1)+[Yy Yx]'*Sw^(-1)*[Yy Yx])^(-1)
%


% Sv = eye(19*NB)./sModel;
% Sy = eye(7*NB)./sUknown;
% Sw = eye(my)./sMeas;
% Dy = Dx\Dy;
% % Sv   = Dx^(-1)*Sv*Dx^(-1)';
% % S    = [Sv + Dy*Sy*Dy', -Dy*Sy; -Sy*Dy', Sy];
% Svinv = Dx'*Sv*Dx;
% Sinv  = [Sv Svinv*Dy; Dy'*Svinv, Sy+ Dy'*Svinv*Dy];
% Sd    = (Sinv+[Yx Yy]'*Sw*[Yx Yy])^(-1);
% Sdx   = Sd(1:NB*19      , 1:NB*19);
% Sdy   = Sd(1+NB*19 : end, 1+NB*19 : end);
% Sdxc = mat2cell(Sdx, 19*ones(1, NB), 19*ones(1, NB));
% Sdyc = mat2cell(Sdy,  7*ones(1, NB),  7*ones(1, NB));

Sv_inv = eye(19*NB)./sModel;
Sy_inv = eye(7*NB) ./sUknown;
Sw_inv = eye(my)   ./sMeas;
% Dy_new = Dx\Dy;
% Sv     = inv(Sv_inv);
% Sy     = inv(Sy_inv);
% Sw     = inv(Sw_inv);

% Sv_new   = Dx^(-1)*Sv*Dx^(-1)';
% S        = [Sv_new + Dy_new*Sy*Dy_new', -Dy_new*Sy; -Sy*Dy_new', Sy];
% Sinv     = S^(-1);
Sinv     = [Dx'*Sv_inv*Dx Dx'*Sv_inv*Dy; Dy'*Sv_inv*Dx, Sy_inv+ Dy'*Sv_inv*Dy];
Sd       = (Sinv+[Yx Yy]'*Sw_inv*[Yx Yy])^(-1);
Sdx      = Sd(1:NB*19      , 1:NB*19);
Sdy      = Sd(1+NB*19 : end, 1+NB*19 : end);
Sdxc     = mat2cell(Sdx, 19*ones(1, NB), 19*ones(1, NB));
Sdyc     = mat2cell(Sdy,  7*ones(1, NB),  7*ones(1, NB));


Sdc  = cell(NB,1);
Sa   = cell(NB,1);
SfB  = cell(NB,1);
Sf   = cell(NB,1);
Stau = cell(NB,1);
Sfx  = cell(NB,1);
Sd2q = cell(NB,1);

for i = 1 : NB
  Sdc{i}   = [Sdxc{i,i}, zeros(19, 7); zeros(7, 19) Sdyc{i,i}];
  Sa{i}    = Sdc{i}(1:6  , 1:6);
  SfB{i}   = Sdc{i}(7:12 , 7:12);
  Sf{i}    = Sdc{i}(13:18, 13:18);
  Stau{i}  = Sdc{i}(19   , 19);
  Sfx{i}   = Sdc{i}(20:25, 20:25);
  Sd2q{i}  = Sdc{i}(26   , 26);
end

