function res = testBNEA(dmodel, ymodel)

res = 0;

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myBNEA    = BNEA(myModel, mySens);
myBNEA    = myBNEA.setEngine('jtree_inf_engine');
disp('[BNEA] JTREE PERFORMANCES...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myPNEA    = PNEA(myModel, mySens);

disp('[BNEA]   - First run -')
tic;
myBNEA = myBNEA.setState(q, dq);
myBNEA = myBNEA.setY(y);
myBNEA = myBNEA.solveID();
t_BNEA = toc;
disp(['[BNEA]   Computation time for BNEA is: ' num2str(t_BNEA) '[sec]']);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[BNEA]   Computation time for PNEA is: ' num2str(t_PNEA) '[sec]']);

disp(['[BNEA]   Diff d  between PNEA and BNEA is ' num2str(norm(myBNEA.d-myPNEA.d))]);
if norm(myBNEA.d-myPNEA.d) > 0.001
   disp('[BNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
disp(['[BNEA]   Diff Sd between PNEA and BNEA is ' num2str(norm(myBNEA.Sd-myPNEA.Sd))]);
disp(['[BNEA]   Diff diag(Sd) between PNEA and BNEA is ' num2str(norm(diag(myBNEA.Sd)-diag(myPNEA.Sd)))]);
if norm(diag(myBNEA.Sd)-diag(myPNEA.Sd)) > 0.1
   disp('[BNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

disp('[BNEA]   - Second run -')
q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myBNEA = myBNEA.setState(q, dq);
myBNEA = myBNEA.setY(y);
myBNEA = myBNEA.solveID();
t_BNEA = toc;
disp(['[BNEA]   Computation time for BNEA is: ' num2str(t_BNEA) '[sec]']);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[BNEA]   Computation time for PNEA is: ' num2str(t_PNEA) '[sec]']);

disp(['[BNEA]   Diff d  between PNEA and BNEA is ' num2str(norm(myBNEA.d-myPNEA.d))]);
if norm(myBNEA.d-myPNEA.d) > 0.001
   disp('[BNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
disp(['[BNEA]   Diff Sd between PNEA and BNEA is ' num2str(norm(myBNEA.Sd-myPNEA.Sd))]);
disp(['[BNEA]   Diff diag(Sd) between PNEA and BNEA is ' num2str(norm(diag(myBNEA.Sd)-diag(myPNEA.Sd)))]);
if norm(diag(myBNEA.Sd)-diag(myPNEA.Sd)) > 0.1
   disp('[BNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

fprintf('\n')
% clear all
% close all

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myBNEA    = BNEA(myModel, mySens);
myBNEA    = myBNEA.setEngine('gaussian_inf_engine');
disp('[BNEA] GAUSS PERFORMANCES...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myPNEA    = PNEA(myModel, mySens);

disp('[BNEA]   - First run -')
tic;
myBNEA = myBNEA.setState(q, dq);
myBNEA = myBNEA.setY(y);
myBNEA = myBNEA.solveID();
t_BNEA = toc;
disp(['[BNEA]   Computation time for BNEA is: ' num2str(t_BNEA) '[sec]']);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[BNEA]   Computation time for PNEA is: ' num2str(t_PNEA) '[sec]']);

disp(['[BNEA]   Diff d  between PNEA and BNEA is ' num2str(norm(myBNEA.d-myPNEA.d))]);
if norm(myBNEA.d-myPNEA.d) > 1e-2
   disp('[BNEA]   Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
disp(['[BNEA]   Diff Sd between PNEA and BNEA is ' num2str(norm(myBNEA.Sd-myPNEA.Sd))]);
if norm(myBNEA.Sd-myPNEA.Sd) > 2
   disp('[BNEA]   Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
disp(['[BNEA]   Diff diag(Sd) between PNEA and BNEA is ' num2str(norm(diag(myBNEA.Sd)-diag(myPNEA.Sd)))]);
if norm(diag(myBNEA.Sd)-diag(myPNEA.Sd)) > 1
   disp('[BNEA]   Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

disp('[BNEA]   - Second run -')
q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myBNEA = myBNEA.setState(q, dq);
myBNEA = myBNEA.setY(y);
myBNEA = myBNEA.solveID();
t_BNEA = toc;
disp(['[BNEA]   Computation time for BNEA is: ' num2str(t_BNEA) '[sec]']);

tic;
myPNEA = myPNEA.setState(q, dq);
myPNEA = myPNEA.setY(y);
myPNEA = myPNEA.solveID();
t_PNEA = toc;
disp(['[BNEA]   Computation time for PNEA is: ' num2str(t_PNEA) '[sec]']);

disp(['[BNEA]   Diff d  between PNEA and BNEA is ' num2str(norm(myBNEA.d-myPNEA.d))]);
if norm(myBNEA.d-myPNEA.d) > 1e-2
   disp('[BNEA]   Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
disp(['[BNEA]   Diff Sd between PNEA and BNEA is ' num2str(norm(myBNEA.Sd-myPNEA.Sd))]);
if norm(myBNEA.Sd-myPNEA.Sd) > 2
   disp('[BNEA]   Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
disp(['[BNEA]   Diff diag(Sd) between PNEA and BNEA is ' num2str(norm(diag(myBNEA.Sd)-diag(myPNEA.Sd)))]);
if norm(diag(myBNEA.Sd)-diag(myPNEA.Sd)) > 1
   disp('[BNEA]   Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

fprintf('\n')
