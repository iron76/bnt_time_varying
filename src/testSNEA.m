clear all
close all
clc

NB        = 5;
m         = 4;
dmodel    = autoTree(NB);
ymodel    = autoSensSNEA(dmodel);

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mySNEA    = SNEA(myModel, mySens);
mySNEA    = mySNEA.setQ(q);
mySNEA    = mySNEA.setDq(dq);
mySNEA    = mySNEA.setY(y);


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
