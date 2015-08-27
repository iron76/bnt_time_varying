function res = testRegressorID

model = autoTree(20,2);
res = 0;
tol = 1e-7;

% create random q,dq,ddq
q = zeros(model.NB,1);
dq = zeros(model.NB,1);
ddq = zeros(model.NB,1);

[tau, a, fB, f] = ID(model, q, dq, ddq);

[fBaseRegr] = regressorID(model, q, dq, ddq);

inertialParams = inertialParametersFromModel(model);

fBaseFromID = f{1};
fBaseFromRegr = fBaseRegr*inertialParams;

if( not(all(isalmost(fBaseFromID,fBaseFromRegr,tol))) )
    disp('Something wrong base force regressors')
    res = 1;
end

end