clear all
close all
clc

NB        = 30;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myPNEA    = PNEA(myModel, mySens);

tic;
myBNEA = myBNEA.setState(q, dq);
myBNEA = myBNEA.setY(y);
myBNEA = myBNEA.solveID();
t_BNEA = toc;
disp(['[1] Computation time for BNEA is: ' num2str(t_BNEA) '[sec]']);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[1] Computation time for PNEA is: ' num2str(t_PNEA) '[sec]']);

disp(['[1] Diff d  between PNEA and BNEA is ' num2str(norm(myBNEA.d-myPNEA.d))]);
disp(['[1] Diff Sd between PNEA and BNEA is ' num2str(norm(myBNEA.Sd-myPNEA.Sd))]);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myBNEA = myBNEA.setState(q, dq);
myBNEA = myBNEA.setY(y);
myBNEA = myBNEA.solveID();
t_BNEA = toc;
disp(['[2] Computation time for BNEA is: ' num2str(t_BNEA) '[sec]']);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[2] Computation time for PNEA is: ' num2str(t_PNEA) '[sec]']);

disp(['[2] Diff d  between PNEA and BNEA is ' num2str(norm(myBNEA.d-myPNEA.d))]);
disp(['[2] Diff Sd between PNEA and BNEA is ' num2str(norm(myBNEA.Sd-myPNEA.Sd))]);

