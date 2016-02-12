function  [tau, a, v,fB, f] = IDv( model, q, qd, qdd, f_ext)

% IDv is a different version of Inverse Dynamics via Recursive Newton-Euler
% Algorithm in Fetherstone toolbox.  It returns also velocities.
% 

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