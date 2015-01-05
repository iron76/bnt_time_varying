clear all
close all
clc

run('iCub.m')
dmodel = iCub_dmodel;
ymodel = autoSensRNEA(dmodel);

dmodel = autoTreeStochastic(dmodel);
ymodel = autoSensStochastic(ymodel);

q      = rand(dmodel.NB,1);
dq     = rand(dmodel.NB,1);
y      = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mySNEA    = SNEA(myModel, mySens);
mySNEA    = mySNEA.setState(q,dq);
mySNEA    = mySNEA.setY(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRNEA    = RNEA(myModel, mySens);
myRNEA    = myRNEA.setState(q,dq);
myRNEA    = myRNEA.setY(y);

if (sum(q-mySNEA.IDstate.q))
   error('Something wrong with the setQ method');
end

if (sum(dq-mySNEA.IDstate.dq))
   error('Something wrong with the setDq method');
end

if (sum(y-mySNEA.IDmeas.y))
   error('Something wrong with the setY method');
end

mySNEA = mySNEA.solveID();
mySNEA.d;

myRNEA = myRNEA.solveID();
myRNEA.d;

disp(['Diff between RNEA and SNEA is ' num2str(norm(mySNEA.d-myRNEA.d))]);

