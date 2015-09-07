function res = testPNEA(dmodel, ymodel)
res = 0;

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

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
disp(['[PNEA] [1st] CPU time for PNEA is: ' num2str(t_PNEA) '[sec]']);

tic;
mySNEA = mySNEA.setState(q,dq);
mySNEA = mySNEA.setY(y);
mySNEA = mySNEA.solveID();
t_SNEA = toc;
disp(['[PNEA] [1st] CPU time for SNEA is: ' num2str(t_SNEA) '[sec]']);
disp(['[PNEA] [1st] Diff SNEA and PNEA is ' num2str(norm(myPNEA.d-mySNEA.d))]);

if norm(myPNEA.d-mySNEA.d) > 1e-9
   disp('[PNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[PNEA] [2st] CPU time for PNEA is: ' num2str(t_PNEA) '[sec]']);

tic;
mySNEA = mySNEA.setState(q,dq);
mySNEA = mySNEA.setY(y);
mySNEA = mySNEA.solveID();
t_SNEA = toc;
disp(['[PNEA] [2st] CPU time for SNEA is: ' num2str(t_SNEA) '[sec]']);
disp(['[PNEA] [2st] Diff SNEA and PNEA is ' num2str(norm(myPNEA.d-mySNEA.d))]);

if norm(myPNEA.d-mySNEA.d) > 1e-9
   disp('[PNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[PNEA] [3st] CPU time for PNEA is: ' num2str(t_PNEA) '[sec]']);

tic;
mySNEA = mySNEA.setState(q,dq);
mySNEA = mySNEA.setY(y);
mySNEA = mySNEA.solveID();
t_SNEA = toc;
disp(['[PNEA] [3st] CPU time for SNEA is: ' num2str(t_SNEA) '[sec]']);
disp(['[PNEA] [3st] Diff SNEA and PNEA is ' num2str(norm(myPNEA.d-mySNEA.d))]);

if norm(myPNEA.d-mySNEA.d) > 1e-9
   disp('[PNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
