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
myMNEA    = MNEA(myModel, mySens);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myPNEA    = PNEA(myModel, mySens);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['Computation time for PNEA is: ' num2str(t_PNEA) '[sec]']);

tic;
myMNEA = myMNEA.setState(q, dq);
myMNEA = myMNEA.setY(y);
myMNEA = myMNEA.solveID();
t_MNEA = toc;
disp(['Computation time for MNEA is: ' num2str(t_MNEA) '[sec]']);

disp(['Diff d  between PNEA and MNEA is ' num2str(norm(myMNEA.d-myPNEA.d))]);
disp(['Diff Sd between PNEA and MNEA is ' num2str(norm(myMNEA.Sd-myPNEA.Sd))]);

