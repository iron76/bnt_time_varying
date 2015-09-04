function res = testMNEA(dmodel, ymodel)

res = 0;

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

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
disp(['[MNEA] Computation time for PNEA is: ' num2str(t_PNEA) '[sec]']);

tic;
myMNEA = myMNEA.setState(q, dq);
myMNEA = myMNEA.setY(y);
myMNEA = myMNEA.solveID();
t_MNEA = toc;
disp(['[MNEA] Computation time for MNEA is: ' num2str(t_MNEA) '[sec]']);

disp(['[MNEA] Diff d  between PNEA and MNEA is ' num2str(norm(myMNEA.d-myPNEA.d))]);
if norm(myMNEA.d-myPNEA.d) > 1e-9
   disp('[MNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
disp(['[MNEA] Diff Sd between PNEA and MNEA is ' num2str(norm(myMNEA.Sd-myPNEA.Sd))]);
if norm(myMNEA.Sd-myPNEA.Sd) ~=0
   disp('[MNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

