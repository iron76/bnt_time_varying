clear all
close all
clc

NB        = 5;
dmodel    = autoTree(NB);

ymodel    = autoSensSNEA(dmodel);

dmodel    = autoTreeStochastic(dmodel);
ymodel    = autoSensStochastic(ymodel);

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myPNEA    = PNEA(myModel, mySens);
myPNEA    = myPNEA.setQ(q);
myPNEA    = myPNEA.setDq(dq);
myPNEA    = myPNEA.setY(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mySNEA    = SNEA(myModel, mySens);
mySNEA    = mySNEA.setQ(q);
mySNEA    = mySNEA.setDq(dq);
mySNEA    = mySNEA.setY(y);

if (sum(q-myPNEA.IDstate.q))
   error('Something wrong with the setQ method');
end

if (sum(dq-myPNEA.IDstate.dq))
   error('Something wrong with the setDq method');
end

if (sum(y-myPNEA.IDmeas.y))
   error('Something wrong with the setY method');
end

myPNEA = myPNEA.solveID();
myPNEA.d;
myPNEA.Sd;

mySNEA = mySNEA.solveID();
mySNEA.d;

disp(['Diff between SNEA and PNEA is ' num2str(norm(myPNEA.d-mySNEA.d))]);

