clear all
close all
clc

run('iCub.m')
dmodel = iCub_dmodel;
ymodel = iCubSens(dmodel);

dmodel = autoTreeStochastic(dmodel);
ymodel = autoSensStochastic(ymodel);

q      = rand(dmodel.NB,1);
dq     = rand(dmodel.NB,1);
y      = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myPNEA    = PNEA(myModel, mySens);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mySNEA    = SNEA(myModel, mySens);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[1st] CPU time for PNEA is: ' num2str(t_PNEA) '[sec]']);

tic;
mySNEA = mySNEA.setState(q,dq);
mySNEA = mySNEA.setY(y);
mySNEA = mySNEA.solveID();
t_SNEA = toc;
disp(['[1st] CPU time for SNEA is: ' num2str(t_SNEA) '[sec]']);
disp(['[1st] Diff SNEA and PNEA is ' num2str(norm(myPNEA.d-mySNEA.d))]);
disp('/\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ ')

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[2st] CPU time for PNEA is: ' num2str(t_PNEA) '[sec]']);

tic;
mySNEA = mySNEA.setState(q,dq);
mySNEA = mySNEA.setY(y);
mySNEA = mySNEA.solveID();
t_SNEA = toc;
disp(['[2st] CPU time for SNEA is: ' num2str(t_SNEA) '[sec]']);
disp(['[2st] Diff SNEA and PNEA is ' num2str(norm(myPNEA.d-mySNEA.d))]);
disp('/\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ ')

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[3st] CPU time for PNEA is: ' num2str(t_PNEA) '[sec]']);

tic;
mySNEA = mySNEA.setState(q,dq);
mySNEA = mySNEA.setY(y);
mySNEA = mySNEA.solveID();
t_SNEA = toc;
disp(['[3st] CPU time for SNEA is: ' num2str(t_SNEA) '[sec]']);
disp(['[3st] Diff SNEA and PNEA is ' num2str(norm(myPNEA.d-mySNEA.d))]);
disp('/\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ ')

