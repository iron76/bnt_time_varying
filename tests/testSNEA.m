function res = testSNEA(dmodel)

res = 0;

ymodel_RNEA  = autoSensRNEA(dmodel);
ymodel_SNEA  = autoSensStochastic( ymodel_RNEA );
mySens_SNEA  = sensors( ymodel_SNEA );

dmodel_SNEA  = autoTreeStochastic(dmodel);
% set Sw_inv to zero
i = dmodel_SNEA.Sw_inv.i;
j = dmodel_SNEA.Sw_inv.j;
for k = 1:length(i)
   dmodel_SNEA.Sw_inv = set(dmodel_SNEA.Sw_inv, zeros(size(dmodel_SNEA.Sw_inv(i(k),j(k)))), i(k), j(k));
end
myModel_SNEA = model(dmodel_SNEA);

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel_RNEA.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRNEA    = RNEA(myModel_SNEA, mySens_SNEA);
myRNEA    = myRNEA.setState(q,dq);
myRNEA    = myRNEA.setY(y);
myRNEA    = myRNEA.solveID();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mySNEA    = SNEA(myModel_SNEA, mySens_SNEA);
mySNEA    = mySNEA.setState(q,dq);
mySNEA    = mySNEA.setY(y);
mySNEA    = mySNEA.solveID();

if (sum(q-mySNEA.IDstate.q))
   disp('[SNEA] Something wrong with the setQ method');
   res = 1;
end

if (sum(dq-mySNEA.IDstate.dq))
   disp('[SNEA] Something wrong with the setDq method');
   res = 1;
end

if (sum(y-mySNEA.IDmeas.y))
   disp('[SNEA] Something wrong with the setY method');
   res = 1;
end


disp(['[SNEA] Diff between RNEA and SNEA is ' num2str(norm(mySNEA.d-myRNEA.d))]);
if norm(mySNEA.d-myRNEA.d) > 1e-3
   disp('[SNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end   

