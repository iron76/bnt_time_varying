clear all
close all
clc

NB        = 20;
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
myBNEA    = BNEA(myModel, mySens);
myBNEA    = myBNEA.setQ(q);
myBNEA    = myBNEA.setDq(dq);
myBNEA    = myBNEA.setY(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myPNEA    = PNEA(myModel, mySens);
myPNEA    = myPNEA.setQ(q);
myPNEA    = myPNEA.setDq(dq);
myPNEA    = myPNEA.setY(y);

if (sum(q-myBNEA.IDstate.q))
   error('Something wrong with the setQ method');
end

if (sum(dq-myBNEA.IDstate.dq))
   error('Something wrong with the setDq method');
end

if (sum(y-myBNEA.IDmeas.y))
   error('Something wrong with the setY method');
end

myBNEA = myBNEA.solveID();
myBNEA.d;
myBNEA.Sd;

myPNEA = myPNEA.solveID();
myPNEA.d;

disp(['Diff d  between PNEA and BNEA is ' num2str(norm(myBNEA.d-myPNEA.d))]);
disp(['Diff Sd between PNEA and BNEA is ' num2str(norm(myBNEA.Sd-myPNEA.Sd))]);

