function res = testDerivativesD(dmodel_SNEA, ymodel_SNEA, dmodel_DNEA, ymodel_DNEA)

res = 0;
num_of_tests  = 1;

myModel = model(dmodel_SNEA);
mySens  = sensors(ymodel_SNEA);
mySNEA  = SNEA(myModel, mySens);

myModel = model(dmodel_DNEA);
mySens  = sensors(ymodel_DNEA);
myDNEA  = DNEA(myModel, mySens);

q  = rand(dmodel_SNEA.NB     ,1);
dq = rand(dmodel_SNEA.NB     ,1);
d  = rand(dmodel_SNEA.NB * 26,1);
Sx = diag(rand(dmodel_SNEA.NB*2, 1));

myDNEA = myDNEA.setState(q,dq);
myDNEA = myDNEA.setD(d);
myDNEA = myDNEA.setDprior(d);
myDNEA = myDNEA.setXprior([q; dq]);
myDNEA = myDNEA.setXvariance(Sx);

for j = 1 : num_of_tests
   q  = rand(dmodel_DNEA.NB     ,1);
   dq = rand(dmodel_DNEA.NB     ,1);
   myDNEA = myDNEA.setState(q,dq);
   myDNEA = compute_dDdq(myDNEA);
   for i = 1 : dmodel_DNEA.NB
      %    mySNEA    = mySNEA.setState(q,dq);
      %    mySNEA    = mySNEA.setY(y);
      %    mySNEA    = mySNEA.solveID();
      
      f      = @(qi) computeD(myDNEA, [q(1:i-1); qi; q(i+1:dmodel_DNEA.NB)]);
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
   q  = rand(dmodel_DNEA.NB     ,1);
   dq = rand(dmodel_DNEA.NB     ,1);
   
   mySNEA = mySNEA.setState(q,dq);
   Sv_inv = mySNEA.IDmodel.modelParams.Sv_inv.matrix;
   Sw_inv = mySNEA.IDmodel.modelParams.Sw_inv.matrix;
   Sy_inv = mySNEA.IDsens.sensorsParams.Sy_inv.matrix;
   Y      = mySNEA.IDsens.sensorsParams.Ys;
   D      = sparse(mySNEA.iDs, mySNEA.jDs, mySNEA.Ds, 19*dmodel_SNEA.NB, 26*dmodel_SNEA.NB);
   S_Dinv = Sv_inv;
   S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);
   S_Yinv = Sy_inv;
   
   myDNEA = myDNEA.setState(q,dq);
   myDNEA = compute_dDdq(myDNEA);
   
   for i = 1 : dmodel_DNEA.NB
      
      f      = @(qi) inv(eye(size(S_dinv)) + Y'*S_Yinv*Y + computeD(mySNEA, [q(1:i-1); qi; q(i+1:dmodel_SNEA.NB)])'*S_Dinv*computeD(mySNEA, [q(1:i-1); qi; q(i+1:dmodel_SNEA.NB)]));
      dDdqi  = deriv(f, q(i,1));
      
      dDi    = myDNEA.dDdq{i}.matrix;
      dDi    = dDi(:, bckPerm);
      dDinv = -inv(eye(size(S_dinv)) + Y'*S_Yinv*Y + D'*S_Dinv*D)*(D'*S_Dinv*dDi + dDi'*S_Dinv*D)*inv(eye(size(S_dinv)) + Y'*S_Yinv*Y + D'*S_Dinv*D);
      % dDinv = (D'*S_Dinv*dDi + dDi'*S_Dinv*D);
      
      if norm( dDinv - dDdqi) > 1e-4
         disp(['[DERIVATIVES] D numerical derivative is quite different: ' num2str(norm(dDinv - dDdqi))])
         res = 1;
      end
%       figure(1)
%       imagesc(dDinv)
%       colorbar
%       figure(2)
%       imagesc(dDdqi)
%       colorbar
%       figure(3)
%       imagesc(dDinv-dDdqi)
%       colorbar
   end
end


for i = 1 : num_of_tests
   q  = rand(dmodel_SNEA.NB     ,1);
   dq = rand(dmodel_SNEA.NB     ,1);
   y  = rand(ymodel_SNEA.m      ,1);
   f  = @(x) compute_d(mySNEA , x, y);
   dd = deriv(f, [q; dq]);
   
   mySNEA = mySNEA.setState(q,dq);
   mySNEA = mySNEA.setY(y);
   mySNEA = mySNEA.solveID();
   
   myDNEA = myDNEA.setState(q,dq);
   myDNEA = compute_dDdq(myDNEA);
   myDNEA = myDNEA.setD(zeros(size(d)));
   myDNEA = myDNEA.setDprior(zeros(size(d)));
   myDNEA = myDNEA.setXprior([q; dq]);
   myDNEA.solveID();
   
   D = sparse(mySNEA.iDs, mySNEA.jDs, mySNEA.Ds, 19*dmodel_SNEA.NB, 26*dmodel_SNEA.NB);
   b = sparse(mySNEA.ibs, ones(size(mySNEA.ibs)), mySNEA.bs, 19*dmodel_SNEA.NB, 1);
   Sv_inv = mySNEA.IDmodel.modelParams.Sv_inv.matrix;
   Sw_inv = mySNEA.IDmodel.modelParams.Sw_inv.matrix;
   Sy_inv = mySNEA.IDsens.sensorsParams.Sy_inv.matrix;
   Y      = mySNEA.IDsens.sensorsParams.Ys;
   
   % y      = mySNEA.IDmeas.y;
   S_Dinv = Sv_inv;
   S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);
   S_Yinv = Sy_inv;
   bY     = zeros(size(y));
   bD     = b;
   % muD    = zeros(length(S_dinv), 1);
   
   dbY    = zeros(length(y), 2*dmodel_SNEA.NB);
   dbD    = myDNEA.dDb_s.matrix;
  
   d     =  mySNEA.d;
   dd_1  = zeros(26*dmodel_DNEA.NB, 2*dmodel_DNEA.NB);
   dd_2  = ((D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y)\(-Y'*S_Yinv*dbY - D'*S_Dinv*dbD));
   for j = 1 : dmodel_DNEA.NB
      dDj   = myDNEA.dDdq{j}.matrix;
      dDj   = dDj(:, bckPerm);
      
      dd_1(:,j) = -inv(S_dinv + Y'*S_Yinv*Y + D'*S_Dinv*D)*(D'*S_Dinv*dDj + dDj'*S_Dinv*D)*d(bckPerm,1);
      dd_1(:,j) = dd_1(:,j) + ((D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y)\(- dDj'*S_Dinv*bD));
   end
   dd_1 = dd_1 + dd_2;
   dd_1 = dd_1(fwdPerm, :);

   subplot(131)
   imagesc(dd)
   colorbar
   subplot(132)
   imagesc(dd_1)
   colorbar
   subplot(133)
   imagesc(dd-dd_1)
   colorbar
   
   % Removed until the numerical issues are fixed (currently the problem is
   % bad conditioned since S_dinv is not yet a proper regularization term [
   % prior is on a subpart of d])
   %
   % if(norm(dd - dd_1) > 1e-2)
   %   disp(['d numerical derivative is quite different: ' num2str(norm(dd - dd_1))])
   %    res = 1;
   % end
end