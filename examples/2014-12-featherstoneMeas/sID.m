function  [tau, a, fB, f, S1] = sID( model, q, qd , Yc , y, S1)

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

nb = 0;
nD = zeros(NB, 1);
for i = 1 : NB
  % b1  = Xup{i}*(-a_grav); 
  % OR
  % b1 = crm(v{i})*vJ;
  nb = nb + 6;
  % b2 = crf(v{i})*model.I{i}*v{i};
  nb = nb + 6;
  
  % Dii = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
  nD(i) = nD(i) + 6 + jn{i}*6;
  % Dii = [model.I{i} -eye(6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})];
  nD(i) = nD(i) + 36 + 6;
  % Dii = [zeros(6,6) eye(6) -eye(6) zeros(6, jn{i}) -inv(Xa{i}') zeros(6, jn{i})];
  nD(i) = nD(i) + 6 + 6 + 36; 
  % Dii = [zeros(jn{i}, 6) zeros(jn{i}, 6) S{i}' -eye(jn{i}) zeros(jn{i}, 6) zeros(jn{i}, jn{i})];
  nD(i) = nD(i) + 6*jn{i} + jn{i};   
  % Dij = [ Xup{i} zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i}) zeros(12+jn{i}, 24+2*jn{i})];
  if model.parent(i) ~= 0
    nD(i) = nD(i) + 36;
  end
  
  ind_j  = find(model.parent == i);
  for j = ind_j
    % Dc{i,j} = [ zeros(12, 24+2*jn{i})
    %     zeros(6,6) zeros(6,6) Xup{j}' zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})
    %     zeros(jn{i}, 24+2*jn{i})];
    nD(i) = nD(i) + 36;
  end
end

pD  = cumsum(nD); 

ib  = zeros(nb, 1);
b1s = zeros(nb, 1);


iD  = zeros(sum(nD), 1);
jD  = zeros(sum(nD), 1);
D1s = zeros(sum(nD), 1);


for i = 1 : NB
  for j = 1 : NB
    inD = [(j-1)*19 (i-1)*19];

    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD11(1:36,i,j) = aa(:);   jD11(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+(7:12));
    aa = aa';   bb = bb';   iD12(1:36,i,j) = aa(:);   jD12(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+(13:18));
    aa = aa';   bb = bb';   iD13(1:36,i,j) = aa(:);   jD13(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+19);
    aa = aa';   bb = bb';   iD14(1: 6,i,j) = aa(:);   jD14(1: 6,i,j) = bb(:);
    
    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD21(1:36,i,j) = aa(:);   jD21(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+(7:12));
    aa = aa';   bb = bb';   iD22(1:36,i,j) = aa(:);   jD22(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+(13:18));
    aa = aa';   bb = bb';   iD23(1:36,i,j) = aa(:);   jD23(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+19);
    aa = aa';   bb = bb';   iD24(1: 6,i,j) = aa(:);   jD24(1: 6,i,j) = bb(:);
    
    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD31(1:36,i,j) = aa(:);   jD31(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+(7:12));
    aa = aa';   bb = bb';   iD32(1:36,i,j) = aa(:);   jD32(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+(13:18));
    aa = aa';   bb = bb';   iD33(1:36,i,j) = aa(:);   jD33(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+19);
    aa = aa';   bb = bb';   iD34(1: 6,i,j) = aa(:);   jD34(1: 6,i,j) = bb(:);
    
    [aa, bb] = meshgrid(inD(1)+19,inD(2)+(1:6));
    aa = aa';   bb = bb';   iD41(1: 6,i,j) = aa(:);   jD41(1: 6,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+19,inD(2)+(7:12));
    aa = aa';   bb = bb';   iD42(1: 6,i,j) = aa(:);   jD42(1: 6,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+19,inD(2)+(13:18));
    aa = aa';   bb = bb';   iD43(1: 6,i,j) = aa(:);   jD43(1: 6,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+19,inD(2)+19);
    aa = aa';   bb = bb';   iD44(1: 1,i,j) = aa;      jD44(1: 1,i,j) = bb(:);

    inD = [(j-1)*19 19*NB+(i-1)*7];
    
    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD15(1:36,i,j) = aa(:);   jD15(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+7);
    aa = aa';   bb = bb';   iD16(1: 6,i,j) = aa(:);   jD16(1: 6,i,j) = bb(:);

    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD25(1:36,i,j) = aa(:);   jD25(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+7);
    aa = aa';   bb = bb';   iD26(1: 6,i,j) = aa(:);   jD26(1: 6,i,j) = bb(:);

    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD35(1:36,i,j) = aa(:);   jD35(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+7);
    aa = aa';   bb = bb';   iD36(1: 6,i,j) = aa(:);   jD36(1: 6,i,j) = bb(:);

    [aa, bb] = meshgrid(inD(1)+19,inD(2)+(1:6));
    aa = aa';   bb = bb';   iD45(1: 6,i,j) = aa(:);   jD45(1: 6,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+19,inD(2)+7);
    aa = aa';   bb = bb';   iD46(1: 1,i,j) = aa(:);   jD46(1: 1,i,j) = bb(:);    
  end
end

pD = [1; pD+1];

for i = 1:model.NB  
  if model.parent(i) == 0
    % a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
    % D1  = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
    % b1  = Xup{i}*(-a_grav);
    D1           = spalloc(6,26,12);
    D1(1:6, 1:6) = -sparse(1:6,1:6,1);
    D1(1:6, 26)  = sparse(S{i});
    b1s((i-1)*12+1: (i-1)*12+6, 1) = Xup{i}*(-a_grav);
    ib((i-1)*12+1: (i-1)*12+6, 1)  = (i-1)*19+1: (i-1)*19+6;

    D1s(pD(i) : pD(i)+5, 1) = -1*ones(6,1);
    iD (pD(i) : pD(i)+5, 1) = iD11(1:6:36,i,i)+[0 1 2 3 4 5]';
    jD (pD(i) : pD(i)+5, 1) = jD11(1:6:36,i,i);
    pD(i) = pD(i) + 6;
    
    D1s(pD(i): pD(i)+jn{i}*6-1, 1) = S{i}(:);
    iD (pD(i): pD(i)+jn{i}*6-1, 1) = iD16(:,i,i);
    jD (pD(i): pD(i)+jn{i}*6-1, 1) = jD16(:,i,i);  
    pD(i) = pD(i) + jn{i}*6;
  else
    % a{i} = ... + S{i}*qdd(i) + crm(v{i})*vJ;
    vJ = S{i}*qd(i);
    % D1 = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
    % b1 = crm(v{i})*vJ;
    D1           = spalloc(6,26,12);
    D1(1:6, 1:6) = -sparse(1:6,1:6,1);
    D1(1:6, 26)  = sparse(S{i});
    
    b1s((i-1)*12+1: (i-1)*12+6, 1)    = crm(v{i})*vJ;
    ib((i-1)*12+1: (i-1)*12+6, 1)     = (i-1)*19+1: (i-1)*19+6;
    
    D1s(pD(i) : pD(i)+5, 1) = -1*ones(6,1);
    iD (pD(i) : pD(i)+5, 1) = iD11(1:6:36,i,i)+[0 1 2 3 4 5]';
    jD (pD(i) : pD(i)+5, 1) = jD11(1:6:36,i,i);
    pD(i) = pD(i) + 6;
    
    D1s(pD(i): pD(i)+jn{i}*6-1, 1) = S{i}(:);
    iD (pD(i): pD(i)+jn{i}*6-1, 1) = iD16(:,i,i);
    jD (pD(i): pD(i)+jn{i}*6-1, 1) = jD16(:,i,i);  
    pD(i) = pD(i) + jn{i}*6;
    % a{i} = Xup{i}*a{model.parent(i)} + ...
    % Dc{i, model.parent(i)} = [ Xup{i} zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})
    %     zeros(12+jn{i}, 24+2*jn{i})];
    
    Dc{i, model.parent(i)} = sparse(i1,j1, Xup{i}(:), 19, 26);
    j = model.parent(i);
    D1s(pD(i): pD(i)+35, 1) = Xup{i}(:);
    iD (pD(i): pD(i)+35, 1) = iD11(:,j,i);
    jD (pD(i): pD(i)+35, 1) = jD11(:,j,i);  
    pD(i) = pD(i) + 36;    

  end
  % fB{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
  % D2 = [model.I{i} -eye(6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})];
  % b2 = crf(v{i})*model.I{i}*v{i};
  
  D2 = sparse([i1' 1:6],[j1' 7:12], [model.I{i}(:); -ones(6,1)], 6, 26);
  b1s((i-1)*12+7: i*12, 1) = crf(v{i})*model.I{i}*v{i};
  ib((i-1)*12+7: i*12, 1)  = (i-1)*19+7: (i-1)*19+12;

  D1s(pD(i) : pD(i)+35, 1) = model.I{i}(:);
  iD (pD(i) : pD(i)+35, 1) = iD21(:,i,i);
  jD (pD(i) : pD(i)+35, 1) = jD21(:,i,i);
  pD(i) = pD(i) + 36;
  
  D1s(pD(i) : pD(i)+5, 1) = -1*ones(6,1);
  iD (pD(i) : pD(i)+5, 1) = iD22(1:6:36,i,i)+[0 1 2 3 4 5]';
  jD (pD(i) : pD(i)+5, 1) = jD22(1:6:36,i,i);
  pD(i) = pD(i) + 6;
  
  Dc{i,i} = [D1; D2];

  % f{i} = fB{i} - Xa{i}' \ f_ext{i};
  % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};
  % D3 = [zeros(6,6) eye(6) -eye(6) zeros(6, jn{i}) -inv(Xa{i}') zeros(6, jn{i})];
  % b3 = zeros(6,1);
  
  A  = -inv(Xa{i}');
  D3 = sparse([1:6 1:6 i5'],[7:12 13:18 j5'], [ones(6,1); -ones(6,1); A(:)], 6, 26);

  D1s(pD(i) : pD(i)+5, 1) = 1*ones(6,1);
  iD (pD(i) : pD(i)+5, 1) = iD32(1:6:36,i,i)+[0 1 2 3 4 5]';
  jD (pD(i) : pD(i)+5, 1) = jD32(1:6:36,i,i);
  pD(i) = pD(i) + 6;

  D1s(pD(i) : pD(i)+5, 1) = -1*ones(6,1);
  iD (pD(i) : pD(i)+5, 1) = iD33(1:6:36,i,i)+[0 1 2 3 4 5]';
  jD (pD(i) : pD(i)+5, 1) = jD33(1:6:36,i,i);
  pD(i) = pD(i) + 6;

  D1s(pD(i) : pD(i)+35, 1) = A(:);
  iD (pD(i) : pD(i)+35, 1) = iD35(:,i,i);
  jD (pD(i) : pD(i)+35, 1) = jD35(:,i,i);
  pD(i) = pD(i) + 36;
    
  % tau(i,1) = S{i}' * f{i};
  % D4 = [zeros(jn{i}, 6) zeros(jn{i}, 6) S{i}' -eye(jn{i}) zeros(jn{i}, 6) zeros(jn{i}, jn{i})];
  % b4 =  zeros(jn{i}, 1);
  
  D4 = sparse([ones(1,6) ones(1,jn{i})], [13:18 19], [S{i}' -ones(1,jn{i})], jn{i}, 26);
  
  D1s(pD(i) : pD(i)+5, 1) = S{i}';
  iD (pD(i) : pD(i)+5, 1) = iD43(:,i,i);
  jD (pD(i) : pD(i)+5, 1) = jD43(:,i,i);
  pD(i) = pD(i) + 6*jn{i};
  
  D1s(pD(i) : pD(i)+jn{i}-1, 1) = -ones(1,jn{i});
  iD (pD(i) : pD(i)+jn{i}-1, 1) = iD44(:,i,i);
  jD (pD(i) : pD(i)+jn{i}-1, 1) = jD44(:,i,i);
  pD(i) = pD(i) + jn{i};

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
    
    D1s(pD(i) : pD(i)+35, 1) = A(:);
    iD (pD(i) : pD(i)+35, 1) = iD33(:,j,i);
    jD (pD(i) : pD(i)+35, 1) = jD33(:,j,i);
    pD(i) = pD(i) + 36;
  end
  
  Dc{i,i} = [Dc{i,i}; D3; D4];
end

yc = y.values;
y  = cell2mat(yc);
b  = sparse(ib, ones(length(ib),1), b1s, 19*NB, 1);
Ds = sparse(iD, jD, D1s, 19*NB, 26*NB); 


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

Dxs = Ds(1:19*NB, 1:19*NB);
Dys = Ds(1:19*NB, 19*NB+1:26*NB);
norm([Dx Dy]-[full(Dxs) full(Dys)])

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
Sv_inv = sparse(1:19*NB, 1:19*NB, 1./sModel);
% Sy_inv = eye(7*NB) ./sUknown;
Sy_inv = sparse(1:7*NB, 1:7*NB, 1./sUknown);
Sy     = sparse(1:7*NB, 1:7*NB, sUknown);
% Sw_inv = eye(my)   ./sMeas;
Sw_inv = sparse(1:my, 1:my, 1./sMeas);
Sinv   = [Dx'*Sv_inv*Dx Dx'*Sv_inv*Dy; Dy'*Sv_inv*Dx, Sy_inv+ Dy'*Sv_inv*Dy];
Dx_inv = Dx\sparse(1:19*NB, 1:19*NB, 1);

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
end
