function res = testDerivativesANEA(dmodel_ANEA, ymodel_ANEA)

res = 0;
num_of_tests  = 1;

NB = dmodel_ANEA.NB;

myModel = model(dmodel_ANEA);
mySens  = sensors(ymodel_ANEA);
myANEA  = ANEA(myModel, mySens);

myModel = model(dmodel_ANEA);
mySens  = sensors(ymodel_ANEA);
myDANEA  = DANEA(myModel, mySens);

q  = rand(NB     ,1);
dq = rand(NB     ,1);
d  = rand(NB * 7,1);
Sx = diag(rand(NB*2, 1));

myDANEA = myDANEA.setState(q,dq);
myDANEA = myDANEA.setD(d);
myDANEA = myDANEA.setDprior(d);
myDANEA = myDANEA.setXprior([q; dq]);
myDANEA = myDANEA.setXvariance(Sx);

for j = 1 : num_of_tests
   q  = rand(NB     ,1);
   dq = rand(NB     ,1);
   myDANEA = myDANEA.setState(q,dq);
   myDANEA = compute_dDdq(myDANEA);
   for i = 1 : NB
      
      f      = @(qi) computeD(myDANEA, [q(1:i-1); qi; q(i+1:NB)]);
      dDdqi  = deriv(f, q(i,1));
      
      if norm(myDANEA.dDdq{i}.matrix - dDdqi) > 1e-7
         disp(['[DERIVATIVES] dD numerical derivative is quite different: ' num2str(norm(myDANEA.dDdq{i}.matrix - dDdqi))])
         res = 1;
      end
      
   end
end

fwdPerm          = myDANEA.id;
bckPerm(fwdPerm) = 1:length(fwdPerm);

for j = 1 : num_of_tests
   q  = rand(NB     ,1);
   dq = rand(NB     ,1);
   
   myANEA = myANEA.setState(q,dq);
   Sv_inv = myANEA.IDmodel.modelParams.Sv_inv.matrix;
   Sw_inv = myANEA.IDmodel.modelParams.Sw_inv.matrix;
   Sy_inv = myANEA.IDsens.sensorsParams.Sy_inv.matrix;
   Y      = myANEA.IDsens.sensorsParams.Ys;
   D      = myANEA.D.matrix;
   S_Dinv = Sv_inv;
   S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);
   S_Yinv = Sy_inv;
   
   myDANEA = myDANEA.setState(q,dq);
   myDANEA = compute_dDdq(myDANEA);
   
   for i = 1 : NB
      
      f      = @(qi) inv(eye(size(S_dinv)) + Y'*S_Yinv*Y + computeD(myANEA, [q(1:i-1); qi; q(i+1:NB)])'*S_Dinv*computeD(myANEA, [q(1:i-1); qi; q(i+1:NB)]));
      dDdqi  = deriv(f, q(i,1));
      
      dDi    = myDANEA.dDdq{i}.matrix;
      dDinv = -inv(eye(size(S_dinv)) + Y'*S_Yinv*Y + D'*S_Dinv*D)*(D'*S_Dinv*dDi + dDi'*S_Dinv*D)*inv(eye(size(S_dinv)) + Y'*S_Yinv*Y + D'*S_Dinv*D);
      
      if norm( dDinv - dDdqi) > 1e-4
         disp(['[DERIVATIVES] D numerical derivative is quite different: ' num2str(norm(dDinv - dDdqi))])
         res = 1;
      end
   end
end


for i = 1 : num_of_tests
   q      = rand(NB     ,1);
   dq     = rand(NB     ,1);
   y      = rand(ymodel_ANEA.m      ,1);
   f      = @(x) compute_d(myANEA , x, y);
   dd_num = deriv(f, [q; dq]);
   
   myANEA = myANEA.setState(q,dq);
   myANEA = myANEA.setY(y);
   myANEA = myANEA.solveID();
   
   myDANEA = myDANEA.setState(q,dq);
   dd      = myDANEA.compute_dq(myANEA.d);

   
   if(norm(dd_num - dd) > 1)
     disp(['[DERIVATIVES] d numerical derivative is quite different: ' num2str(norm(dd_num - dd))])
      res = 1;
   end
end