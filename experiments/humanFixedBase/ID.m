function  [tau, a, fB, f, v,fReactionBase] = ID( model, q, qd, qdd, f_ext)

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
    f{model.parent(i)} = fB{i} + Xup{i}'*f{i};
  end
end
   % fBase =  zeros(size(f{1}));%model.FootI *(-a_grav) + 
    fReactionBase = f{1} ;%+ Xup{1}'*f{1} ;
    fReactionBase;
end