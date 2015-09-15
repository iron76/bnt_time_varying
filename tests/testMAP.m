function res = testMAP(dmodel, ymodel)
res = 0;

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
d         = rand(26*dmodel.NB,1);
y         = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myMAP    = MAP(myModel, mySens);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mydMAP   = dMAP(myModel, mySens);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mySNEA    = SNEA(myModel, mySens);

tic;
myMAP = myMAP.setState(q, dq);
myMAP = myMAP.setY(y);
myMAP = myMAP.solveID();
t_MAP = toc;
disp(['[MAP]  [1st] CPU time for  MAP is: ' num2str(t_MAP) '[sec]']);

tic;
mydMAP = mydMAP.setState(q, dq);
mydMAP = mydMAP.setY(y);
mydMAP = mydMAP.solveID();
mydMAP = mydMAP.updatedDbSubmatrix(mydMAP.d);
dd_dq  = mydMAP.compute_dq(mydMAP.d);
t_dMAP = toc;
disp(['[dMAP] [1st] CPU time for dMAP is: ' num2str(t_dMAP) '[sec]']);


tic;
mySNEA = mySNEA.setState(q,dq);
mySNEA = mySNEA.setY(y);
mySNEA = mySNEA.solveID();
t_SNEA = toc;
disp(['[MAP]  [1st] CPU time for SNEA is: ' num2str(t_SNEA) '[sec]']);

if norm(myMAP.d-mySNEA.d) ~=0
   disp('[MAP] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

if norm(mydMAP.d-mySNEA.d) ~=0
   disp('[dMAP] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end


q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
d         = rand(26*dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myMAP = myMAP.setState(q, dq);
myMAP = myMAP.setY(y);
myMAP = myMAP.solveID();
t_MAP = toc;
disp(['[MAP]  [2nd] CPU time for  MAP is: ' num2str(t_MAP) '[sec]']);

tic;
mydMAP = mydMAP.setState(q, dq);
mydMAP = mydMAP.setY(y);
mydMAP = mydMAP.solveID();
mydMAP = mydMAP.updatedDbSubmatrix(mydMAP.d);
dd_dq  = mydMAP.compute_dq(mydMAP.d);
t_dMAP = toc;
disp(['[dMAP] [2nd] CPU time for dMAP is: ' num2str(t_dMAP) '[sec]']);


tic;
mySNEA = mySNEA.setState(q,dq);
mySNEA = mySNEA.setY(y);
mySNEA = mySNEA.solveID();
t_SNEA = toc;
disp(['[MAP]  [2nd] CPU time for SNEA is: ' num2str(t_SNEA) '[sec]']);

if norm(myMAP.d-mySNEA.d) ~=0
   disp('[MAP] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

if norm(mydMAP.d-mySNEA.d) ~=0
   disp('[dMAP] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
d         = rand(26*dmodel.NB,1);
y         = rand(ymodel.m,1);

tic;
myMAP = myMAP.setState(q, dq);
myMAP = myMAP.setY(y);
myMAP = myMAP.solveID();
t_MAP = toc;
disp(['[MAP]  [3rd] CPU time for  MAP is: ' num2str(t_MAP) '[sec]']);

tic;
mydMAP = mydMAP.setState(q, dq);
mydMAP = mydMAP.setY(y);
mydMAP = mydMAP.solveID();
mydMAP = mydMAP.updatedDbSubmatrix(mydMAP.d);
dd_dq  = mydMAP.compute_dq(mydMAP.d);
t_dMAP = toc;
disp(['[dMAP] [3rd] CPU time for dMAP is: ' num2str(t_dMAP) '[sec]']);


tic;
mySNEA = mySNEA.setState(q,dq);
mySNEA = mySNEA.setY(y);
mySNEA = mySNEA.solveID();
t_SNEA = toc;
disp(['[MAP]  [3rd] CPU time for SNEA is: ' num2str(t_SNEA) '[sec]']);

if norm(myMAP.d-mySNEA.d) ~=0
   disp('[MAP] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end

if norm(mydMAP.d-mySNEA.d) ~=0
   disp('[dMAP] Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end