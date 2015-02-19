clear all
close all
clc

NB        = 20;
dmodel    = autoTree(NB);
ymodel    = autoSensRNEA(dmodel);

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

if (sum(q-myRNEA.IDstate.q) ~= 0);
   error('Something wrong with the setQ method');
end

if (sum(dq-myRNEA.IDstate.dq) ~= 0);
   error('Something wrong with the setDq method');
end

if (sum(y-myRNEA.IDmeas.y) ~= 0);
   error('Something wrong with the setY method');
end

myRNEA    = myRNEA.solveID();
d_myRNEA  = myRNEA.d;

d2q   = zeros(NB,1);
fx    = cell(1,NB);
for i = 1 : NB
   fx{i}    = y((1:6)+(i-1)*7,1);
   d2q(i,1) = y(7*i); 
end

[tau, a, fB, f] = ID( dmodel, q, dq, d2q, fx);

d = zeros(26*NB, 1);
for i = 1 : NB
   d((1:26)+(i-1)*26, 1) = [a{i}; fB{i}; f{i}; tau(i,1); fx{i}; d2q(i,1)];
end

if (sum(d-d_myRNEA) ~= 0)
   error('Something wrong with the solveID method');
end


[a_myRNEA,fB_myRNEA,f_myRNEA,tau_myRNEA,fx_myRNEA,d2q_myRNEA] = extractDynVar(NB, d_myRNEA);

for i = 1:NB    
    if (sum(a{1,i}-a_myRNEA{1,i}) ~= 0)
        error('Something wrong with the a variable check');
    end
end

for i = 1:NB 
    if (sum(fB{1,i}-fB_myRNEA{1,i}) ~= 0)
        error('Something wrong with the fB variable check');
    end
end

for i = 1:NB 
    if (sum(f{1,i}-f_myRNEA{1,i}) ~= 0)
        error('Something wrong with the f variable check');
    end
end

if (sum(tau-tau_myRNEA) ~= 0)
    error('Something wrong with the tau variable check');
end

for i = 1:NB 
    if (sum(fx{1,i}-fx_myRNEA{1,i}) ~= 0)
        error('Something wrong with the f variable check');
    end
end

if (sum(d2q-d2q_myRNEA) ~= 0)
    error('Something wrong with the f variable check');
end

