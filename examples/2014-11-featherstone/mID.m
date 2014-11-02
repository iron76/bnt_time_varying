function  [tau, a, fB, f] = mID( model, q, qd, qdd, f_ext )

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

sKn = 1e-5;
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

dyc = cell(NB,1);
Dxc = cell(NB,NB);
Dyc = cell(NB,NB);
for i = 1 : NB
  dyc{i,1} = [f_ext{i}; qdd(i)];
  for j = 1 : NB
    Dxc{i,j} = Dc{i,j}(:, 1:19);
    Dyc{i,j} = Dc{i,j}(:, 20:end);
  end
end
dy = cell2mat(dyc);
Dx = cell2mat(Dxc);
Dy = cell2mat(Dyc);
b  = cell2mat(bc);
dx = -inv(Dx)*(b + Dy*dy);
% dx = -Dx\(b + Dy*dy);


dxc   = mat2cell(dx, ones(NB,1).*19, 1);
a   = zeros(NB,6);
fB  = zeros(NB,6);
f   = zeros(NB,6);
tau = zeros(NB,1);

for i = 1 : NB
  a( i,1:6) = dxc{i,1}( 1: 6);
  fB(i,1:6) = dxc{i,1}( 7: 12);
  f( i,1:6) = dxc{i,1}(13: 18);  
  tau(i,1)  = dxc{i,1}(19);
end
