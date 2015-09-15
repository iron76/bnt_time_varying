function res = testConditioning(dmodel, ymodel)

res = 0;

q      = rand(dmodel.NB,1);
dq     = rand(dmodel.NB,1);
y      = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mySNEA = SNEA(myModel, mySens);
mySNEA = mySNEA.setState(q,dq);

Sv_inv = mySNEA.IDmodel.modelParams.Sv_inv.matrix;
Sw_inv = mySNEA.IDmodel.modelParams.Sw_inv.matrix;
Sy_inv = mySNEA.IDsens.sensorsParams.Sy_inv.matrix;
Y      = mySNEA.IDsens.sensorsParams.Ys;
D      = sparse(mySNEA.iDs, mySNEA.jDs, mySNEA.Ds, 19*dmodel.NB, 26*dmodel.NB);
S_Dinv = Sv_inv;
S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);
S_Yinv = Sy_inv;

A = full(D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y);
% A = full(D'*D + Y'*Y);
disp(['[CONDITIONING] cond(A) is ' num2str(cond(A), '%10.5e\n')]);
disp(['[CONDITIONING] Precision in the computations will be: 1e', num2str(round(log10(eps)+log10(cond(A))))])
B = full(diag(1./max(abs(A')))*A);
disp(['Diff between cond(A) and cond(B) is ' num2str(cond(A), '%10.5e\n') ',' num2str(cond(B), '%10.5e\n')]);
% if norm(mySNEA.d-myRNEA.d) > 2.0
%    disp('Result is excessively inaccurate. Test is declared failed!');
%    res = 1;
% end

