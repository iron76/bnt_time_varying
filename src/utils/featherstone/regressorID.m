function  [ fBaseRegressor ] = regressorID( model, q, qd, qdd)

% regressorID  Inverse Dynamics in regressor form
% regressorID(model,q,qd,qdd) calculates the inverse dynamics of a kinematic
% tree in regressor form.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables
% this class is in implemented in a way such that :
%  [tau, a, fB, f] = ID( model, q, qd, qdd)
%  [fBaseRegr, tauRegr] = regressorID(model,q,qd,qdd)
%  inertialParameters = inertialParametersFromModel(model);
% tauRegr*inertialParameters == tau
% fBaseRegressor*inertialParameters == f{1}
% TODO : implement torqueRegressor

a_grav = get_gravity(model);

%tauRegressor    = zeros(model.NB,10*model.NB);
fBaseRegressor  = zeros(6,10*model.NB);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  [~, jn{i}] = size(S{i});
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
    XwrtLinkOne{i} = eye(6,6);
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
    XwrtLinkOne{i} = Xup{i}*XwrtLinkOne{model.parent(i)};
  end
  %fB{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
  fBaseRegressor(:,(1+(i-1)*10):(i*10)) = XwrtLinkOne{i}'*netWrenchRegressor(a{i},v{i});
end