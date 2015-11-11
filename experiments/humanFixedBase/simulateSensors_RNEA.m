function  [ y_RNEA ] = simulateSensors_RNEA( model, q, qd, qdd, f_ext)

% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  [~, jn{i}] = size(S{i});
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
  end
  fB{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
end

if nargin >= 5
  [f, Xa] = apply_external_forces( model.parent, Xup, fB, f_ext );
else
    f = fB;
end


for i = model.NB:-1:1
  tau(i,1) = S{i}' * f{i};
  if model.parent(i) ~= 0
    f{model.parent(i)} = f{model.parent(i)} + Xup{i}'*f{i};
  end
end

a_2_2 = a{2}(4:6)';
%v_2_2 = v{2}(1:3)';
fts(4:6) = f{1}(1:3); 
fts(1:3) = f{1}(4:6);

y_RNEA = zeros(26,1);
y_RNEA(1:6) = fts;
y_RNEA(7:9) = a_2_2;
y_RNEA(10:12) = zeros(3,1);
y_RNEA(13:24) = zeros(12,1);
y_RNEA(25:26) = qdd;