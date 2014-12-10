clear all
close all
clc

NB        = 5;
m         = 4;
dmodel    = autoTree(NB);
ymodel    = autoSensSNEA(NB);

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

for i = 1 : NB
   fx{i}    = y((1:6)+(i-1)*7,1);
   d2q(i,1) = y(7*i);
end

[tau, a, fB, f] = ID( dmodel, q, dq, d2q, fx);

d = zeros(26*NB, 1);
for i = 1 : NB
   d((1:26)+(i-1)*26, 1) = [a{i}; fB{i}; f{i}; tau(i,1); fx{i}; d2q(i,1)];
end

if (sum(d-mySNEA.d))
   error('Something wrong with the solveID method');
end