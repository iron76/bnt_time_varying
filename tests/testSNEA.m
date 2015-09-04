function res = testSNEA(dmodel, ymodel)

res = 0;

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRNEA    = RNEA(myModel, mySens);
myRNEA    = myRNEA.setState(q,dq);
myRNEA    = myRNEA.setY(y);
myRNEA    = myRNEA.solveID();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mySNEA    = SNEA(myModel, mySens);
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
if norm(mySNEA.d-myRNEA.d) > 2.0
   disp('[SNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end   

