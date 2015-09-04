function res = testMRNEA(dmodel, ymodel)

res = 0;

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
disp(['[MRNEA] [1st] CPU time for MRNEA is: ' num2str(t_MRNEA) '[sec]']);

tic;
myRNEA = myRNEA.setState(q,dq);
myRNEA = myRNEA.setY(y);
myRNEA = myRNEA.solveID();
t_RNEA = toc;
disp(['[MRNEA] [1st] CPU time for RNEA is: ' num2str(t_RNEA) '[sec]']);
disp(['[MRNEA] [1st] Diff RNEA and MRNEA is ' num2str(norm(myMRNEA.d-myRNEA.d))]);

if norm(norm(myMRNEA.d-myRNEA.d)) > 1e-10
   disp('[MRNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myMRNEA = myMRNEA.setState(q, dq);
myMRNEA = myMRNEA.setY(y);
myMRNEA = myMRNEA.solveID();
t_MRNEA = toc;
disp(['[MRNEA] [2st] CPU time for MRNEA is: ' num2str(t_MRNEA) '[sec]']);

tic;
myRNEA = myRNEA.setState(q,dq);
myRNEA = myRNEA.setY(y);
myRNEA = myRNEA.solveID();
t_RNEA = toc;
disp(['[MRNEA] [2st] CPU time for RNEA is: ' num2str(t_RNEA) '[sec]']);
disp(['[MRNEA] [2st] Diff RNEA and MRNEA is ' num2str(norm(myMRNEA.d-myRNEA.d))]);

if norm(norm(myMRNEA.d-myRNEA.d)) > 1e-10
   disp('[MRNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myMRNEA = myMRNEA.setState(q, dq);
myMRNEA = myMRNEA.setY(y);
myMRNEA = myMRNEA.solveID();
t_MRNEA = toc;
disp(['[MRNEA] [3st] CPU time for MRNEA is: ' num2str(t_MRNEA) '[sec]']);

tic;
myRNEA = myRNEA.setState(q,dq);
myRNEA = myRNEA.setY(y);
myRNEA = myRNEA.solveID();
t_RNEA = toc;
disp(['[MRNEA] [3st] CPU time for RNEA is: ' num2str(t_RNEA) '[sec]']);
disp(['[MRNEA] [3st] Diff RNEA and MRNEA is ' num2str(norm(myMRNEA.d-myRNEA.d))]);

if norm(norm(myMRNEA.d-myRNEA.d)) > 1e-10
   disp('[MRNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
