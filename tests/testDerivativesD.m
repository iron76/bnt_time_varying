function res = testDerivativesD(dmodel_SNEA, ymodel_SNEA)

res = 0;
num_of_tests  = 1;

NB = dmodel_SNEA.NB;

myModel = model(dmodel_SNEA);
mySens  = sensors(ymodel_SNEA);
mySNEA  = SNEA(myModel, mySens);

myModel = model(dmodel_SNEA);
mySens  = sensors(ymodel_SNEA);
myDNEA  = DNEA(myModel, mySens);

q  = rand(NB     ,1);
dq = rand(NB     ,1);
d  = rand(NB * 26,1);
Sx = diag(rand(NB*2, 1));

myDNEA = myDNEA.setState(q,dq);
myDNEA = myDNEA.setD(d);
myDNEA = myDNEA.setDprior(d);
myDNEA = myDNEA.setXprior([q; dq]);
myDNEA = myDNEA.setXvariance(Sx);

for j = 1 : num_of_tests
   q  = rand(NB     ,1);
   dq = rand(NB     ,1);
   myDNEA = myDNEA.setState(q,dq);
   myDNEA = compute_dDdq(myDNEA);
   for i = 1 : NB
      
      f      = @(qi) computeD(myDNEA, [q(1:i-1); qi; q(i+1:NB)]);
      dDdqi  = deriv(f, q(i,1));
      
      if norm(myDNEA.dDdq{i}.matrix - dDdqi(:, myDNEA.id)) > 1e-7
         disp(['[DERIVATIVES] dD numerical derivative is quite different: ' num2str(norm(myDNEA.dDdq{i}.matrix - dDdqi))])
         res = 1;
      end
      
   end
end

fwdPerm          = myDNEA.id;
bckPerm(fwdPerm) = 1:length(fwdPerm);

for j = 1 : num_of_tests
   q  = rand(NB     ,1);
   dq = rand(NB     ,1);
   
   mySNEA = mySNEA.setState(q,dq);
   Sv_inv = mySNEA.IDmodel.modelParams.Sv_inv.matrix;
   Sw_inv = mySNEA.IDmodel.modelParams.Sw_inv.matrix;
   Sy_inv = mySNEA.IDsens.sensorsParams.Sy_inv.matrix;
   Y      = mySNEA.IDsens.sensorsParams.Ys;
   D      = sparse(mySNEA.iDs, mySNEA.jDs, mySNEA.Ds, 19*NB, 26*NB);
   S_Dinv = Sv_inv;
   S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);
   S_Yinv = Sy_inv;
   
   myDNEA = myDNEA.setState(q,dq);
   myDNEA = compute_dDdq(myDNEA);
   
   for i = 1 : NB
      
      f      = @(qi) inv(eye(size(S_dinv)) + Y'*S_Yinv*Y + computeD(mySNEA, [q(1:i-1); qi; q(i+1:NB)])'*S_Dinv*computeD(mySNEA, [q(1:i-1); qi; q(i+1:NB)]));
      dDdqi  = deriv(f, q(i,1));
      
      dDi    = myDNEA.dDdq{i}.matrix;
      dDi    = dDi(:, bckPerm);
      dDinv = -inv(eye(size(S_dinv)) + Y'*S_Yinv*Y + D'*S_Dinv*D)*(D'*S_Dinv*dDi + dDi'*S_Dinv*D)*inv(eye(size(S_dinv)) + Y'*S_Yinv*Y + D'*S_Dinv*D);
      
      if norm( dDinv - dDdqi) > 1e-4
         disp(['[DERIVATIVES] D numerical derivative is quite different: ' num2str(norm(dDinv - dDdqi))])
         res = 1;
      end
   end
end


for i = 1 : num_of_tests
   q  = rand(NB     ,1);
   dq = rand(NB     ,1);
   y  = rand(ymodel_SNEA.m      ,1);
   f  = @(x) compute_d(mySNEA , x, y);
   dd = deriv(f, [q; dq]);
   
   mySNEA = mySNEA.setState(q,dq);
   mySNEA = mySNEA.setY(y);
   mySNEA = mySNEA.solveID();
   
   myDNEA = myDNEA.setState(q,dq);
   dd_dq  = myDNEA.compute_dq(mySNEA.d);
   
   
   if(norm(dd - dd_dq) > 1e-2)
      disp(['[DERIVATIVES] d numerical derivative is quite different: ' num2str(norm(dd - dd_dq))])
      
      subplot(131)
      imagesc(dd)
      colorbar
      subplot(132)
      imagesc(dd_dq)
      colorbar
      subplot(133)
      imagesc(dd-dd_dq)
      colorbar
      res = 1;
   end
end