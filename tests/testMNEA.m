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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myMAP     =  MAP(myModel, mySens);

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

tic;
myMAP = myMAP.setState(q, dq);
myMAP = myMAP.setY(y);
myMAP = myMAP.solveID('variance');
t_MAP = toc;
disp(['[myMAP]  Computation time for MNEA is: ' num2str(t_MAP) '[sec]']);


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


disp(['[MNEA] Diff d  between MAP  and MNEA is ' num2str(norm(myMAP.d-myMNEA.d))]);
if norm(myMAP.d-myMNEA.d) > 1e-9
   disp('[MNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
disp(['[MNEA] Diff Sd between MAP  and MNEA is ' num2str(norm(myMAP.Sd-myMNEA.Sd))]);
if norm(myMAP.Sd-myMNEA.Sd) ~=0
   disp('[MNEA] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end
