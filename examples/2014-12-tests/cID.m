function  [d, Sd, S1] = cID( model, q, qd , Yc , y, S1)

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

% indexes for writing sparse matrix data
[i1, j1] = meshgrid(1:6,1:6);
i1 = i1';     j1 = j1';
i1 = i1(:);   j1 = j1(:);

[i2, j2] = meshgrid(1:6,7:12);
i2 = i2';     j2 = j2';
i2 = i2(:);   j2 = j2(:);

[i3, j3] = meshgrid(1:6,13:18);
i3 = i3';     j3 = j3';
i3 = i3(:);   j3 = j3(:);

[i5, j5] = meshgrid(1:6,20:25);
i5 = i5';     j5 = j5';
i5 = i5(:);   j5 = j5(:);


for i = 1:model.NB
  if model.parent(i) == 0
    % a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
    % D1  = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
    % b1  = Xup{i}*(-a_grav);
    D1           = spalloc(6,26,12);
    D1(1:6, 1:6) = -sparse(1:6,1:6,1);
    D1(1:6, 26)  = sparse(S{i});
    b1           = sparse(1:6,ones(1,6), Xup{i}*(-a_grav), 6, 1);
  else
    % a{i} = ... + S{i}*qdd(i) + crm(v{i})*vJ;
    vJ = S{i}*qd(i);
    % D1 = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
    % b1 = crm(v{i})*vJ;
    D1           = spalloc(6,26,12);
    D1(1:6, 1:6) = -sparse(1:6,1:6,1);
    D1(1:6, 26)  = sparse(S{i});
    b1           = sparse(1:6,ones(1,6), crm(v{i})*vJ, 6, 1);
    % a{i} = Xup{i}*a{model.parent(i)} + ...
    % Dc{i, model.parent(i)} = [ Xup{i} zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})
    %     zeros(12+jn{i}, 24+2*jn{i})];
    
    Dc{i, model.parent(i)} = sparse(i1,j1, Xup{i}(:), 19, 26);
  end
  % fB{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
  % D2 = [model.I{i} -eye(6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})];
  % b2 = crf(v{i})*model.I{i}*v{i};
  
  D2 = sparse([i1' 1:6],[j1' 7:12], [model.I{i}(:); -ones(6,1)], 6, 26);
  b2 = sparse(1:6,ones(1,6), crf(v{i})*model.I{i}*v{i}         , 6, 1);
  
  Dc{i,i} = [D1; D2];
  bc{i}   = [b1; b2];
end

for i = model.NB:-1:1
  % f{i} = fB{i} - Xa{i}' \ f_ext{i};
  % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};
  % D3 = [zeros(6,6) eye(6) -eye(6) zeros(6, jn{i}) -inv(Xa{i}') zeros(6, jn{i})];
  % b3 = zeros(6,1);
  
  A  = -inv(Xa{i}');
  D3 = sparse([1:6 1:6 i5'],[7:12 13:18 j5'], [ones(6,1); -ones(6,1); A(:)], 6, 26);
  b3 = sparse([],[],[],6,1,0);
  
  ind_j  = find(model.parent == i);
  for j = ind_j
    % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};
    % Dc{i,j} = [ zeros(12, 24+2*jn{i})
    %     zeros(6,6) zeros(6,6) Xup{j}' zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})
    %     zeros(jn{i}, 24+2*jn{i})];
    
    [ii, jj] = meshgrid(13:18,13:18);
    ii = ii';     jj = jj';
    ii = ii(:);   jj = jj(:);
    A       = Xup{j}';
    Dc{i,j} = sparse(ii, jj, A(:), 18+jn{i}, 24+2*jn{i});
  end
  
  % tau(i,1) = S{i}' * f{i};
  % D4 = [zeros(jn{i}, 6) zeros(jn{i}, 6) S{i}' -eye(jn{i}) zeros(jn{i}, 6) zeros(jn{i}, jn{i})];
  % b4 =  zeros(jn{i}, 1);
  
  D4 = sparse([ones(1,6) ones(1,jn{i})], [13:18 19], [S{i}' -ones(1,jn{i})], jn{i}, 26);
  b4 = sparse([],[],[],jn{i},1,0);
  
  Dc{i,i} = [Dc{i,i}; D3; D4];
  bc{i}   = [bc{i};   b3; b4];
end

yc = y.values;
y  = cell2mat(yc);
b  = cell2mat(bc);

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

if nargin <= 5
  Ss = Sinv+[Yx Yy]'*Sw_inv*[Yx Yy];
  [L,~,S1] = chol(Ss, 'lower');    % S1'*W*S1 = L*L'
  % PWinv1 = S1*inv_chol(L)*S1';
  % Ls = S1*L;
else
  Ss = Sinv+[Yx Yy]'*Sw_inv*[Yx Yy];
  L = chol(S1'*Ss*S1, 'lower');    % S1'*W*S1 = L*L'
  % PWinv1 = S1*inv_chol(L)*S1';
  % Ls = S1*L;
end

% Sd = (Ls'\(Ls\eye(NB*26)));
Sd = S1*(L'\(L\eye(NB*26)))*S1';

Sdx   = Sd(1:NB*19      , 1:NB*19);
Sdy   = Sd(1+NB*19 : end, 1+NB*19 : end);


Sdxc = mat2cell(Sdx, 19*ones(1, NB), 19*ones(1, NB));
Sdyc = mat2cell(Sdy,  7*ones(1, NB),  7*ones(1, NB));

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

%%
%
% Solution computation
Dx_inv = inv(Dx);
Sy     = inv(Sy_inv);
Sxy = [Dx_inv + Dx_inv*Dy*Sy*Dy'*Sv_inv, -Dx_inv*Dy*Sy; -Sy*Dy'*Sv_inv, Sy];
mx  = -b;
my  = zeros(7*NB,1);
mxy = -Sxy*[-mx; -Dy'*Sv_inv*mx - Sy_inv*my];
% d   = mxy + Sd*[Yx Yy]'*Sw_inv*(y-[Yx Yy]*mxy);
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

