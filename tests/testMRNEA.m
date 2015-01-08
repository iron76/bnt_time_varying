clear all
close all
clc

NB        = 20;
m         = 4;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myMRNEA   = MRNEA(myModel, mySens);

tic;
myMRNEA = myMRNEA.setState(q, dq);
myMRNEA = myMRNEA.setY(y);
myMRNEA = myMRNEA.solveID();
t_MRNEA = toc;
disp(['[1st] CPU time for MRNEA is: ' num2str(t_MRNEA) '[sec]']);

tic;
myRNEA = myRNEA.setState(q,dq);
myRNEA = myRNEA.setY(y);
myRNEA = myRNEA.solveID();
t_RNEA = toc;
disp(['[1st] CPU time for RNEA is: ' num2str(t_RNEA) '[sec]']);
disp(['[1st] Diff RNEA and MRNEA is ' num2str(norm(myMRNEA.d-myRNEA.d))]);
disp('/\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ ')

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myMRNEA = myMRNEA.setState(q, dq);
myMRNEA = myMRNEA.setY(y);
myMRNEA = myMRNEA.solveID();
t_MRNEA = toc;
disp(['[2st] CPU time for MRNEA is: ' num2str(t_MRNEA) '[sec]']);

tic;
myRNEA = myRNEA.setState(q,dq);
myRNEA = myRNEA.setY(y);
myRNEA = myRNEA.solveID();
t_RNEA = toc;
disp(['[2st] CPU time for RNEA is: ' num2str(t_RNEA) '[sec]']);
disp(['[2st] Diff RNEA and MRNEA is ' num2str(norm(myMRNEA.d-myRNEA.d))]);
disp('/\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ ')

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myMRNEA = myMRNEA.setState(q, dq);
myMRNEA = myMRNEA.setY(y);
myMRNEA = myMRNEA.solveID();
t_MRNEA = toc;
disp(['[3st] CPU time for MRNEA is: ' num2str(t_MRNEA) '[sec]']);

tic;
myRNEA = myRNEA.setState(q,dq);
myRNEA = myRNEA.setY(y);
myRNEA = myRNEA.solveID();
t_RNEA = toc;
disp(['[3st] CPU time for RNEA is: ' num2str(t_RNEA) '[sec]']);
disp(['[3st] Diff RNEA and MRNEA is ' num2str(norm(myMRNEA.d-myRNEA.d))]);
disp('/\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ ')

