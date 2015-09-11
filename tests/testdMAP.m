function res = testdMAP(dmodel, ymodel)

res = 0;

myModel = model(dmodel);
mySens  = sensors(ymodel);
myMAP  = dMAP(myModel, mySens);


q  = rand(dmodel.NB     ,1);
dq = rand(dmodel.NB     ,1);
d  = rand(dmodel.NB * 26,1);
y  = rand(ymodel.m      ,1);

myMAP = myMAP.setState(q, dq);
myMAP = myMAP.setY(y);
myMAP = myMAP.solveID();
myMAP = myMAP.updatedDbSubmatrix(myMAP.d);
dd_dq  = myMAP.compute_dq(myMAP.d);

f  = @(x) computeDb(myMAP , myMAP.d , x);
dD = deriv(f, [q; dq]);

f  = @(x) compute_d(myMAP , x, y);
dd = deriv(f, [q; dq]);


if norm(myMAP.dDd_s.matrix+myMAP.dbD_s.matrix - dD) > 1e-3
   disp(['[DERIVATIVES] dD numerical derivative is quite different: ' num2str(norm(myMAP.dDd_s.matrix+myMAP.dbD_s.matrix - dD))])
   res = 1;
end

subplot(131)
imagesc(dd)
colorbar
subplot(132)
imagesc(dd_dq)
colorbar
subplot(133)
imagesc(dd-dd_dq)
colorbar




