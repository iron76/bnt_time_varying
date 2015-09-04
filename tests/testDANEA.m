function res = testDANEA(dmodel, ymodel)

res = 0;

NB            = dmodel.NB;
num_of_tests  = 1;
ymodel_DANEA  = autoSensDANEA(dmodel, ymodel, zeros(NB,1), zeros(NB, 1));
myANEA        = DANEA(model(dmodel), sensors(ymodel_DANEA));


for i = 1 : num_of_tests
   q  = rand(dmodel.NB     ,1);
   dq = rand(dmodel.NB     ,1);
   d  = rand(dmodel.NB * 7,1);
   
   myANEA = myANEA.setState(q,dq);
   
   f  = @(x) computeY(myANEA, d, x);
   dy = deriv(f, [q; dq]);
   
   % In current formulation 
   %
   % y = Yd d + Yx x + b_Y(x) = [Yd(x) Yx][d ; x] + b_Y(x)
   %
   % and therefore the derivative is:
   %
   % dy/dx = Yx + db_Y/dx
   
   Yx  = myANEA.IDsens.sensorsParams.Ys(:,(end-2*dmodel.NB+1):end);
   dby = myANEA.dby_s.matrix + Yx;
   if norm( dby - dy) > 1e-6
      disp(['[DANEA] dy numerical derivative is quite different: ' num2str(norm(dby - dy))])
      imagesc([dby, dy])
      colorbar
      res = 1;
   end
end


for i = 1 : num_of_tests
   q  = rand(dmodel.NB     ,1);
   dq = rand(dmodel.NB     ,1);
   d  = rand(dmodel.NB * 7,1);
   Sx = diag(rand(dmodel.NB*2, 1));
   f  = @(x) computeDb(myANEA , d , x);
   dD = deriv(f, [q; dq]);
   
   myANEA = myANEA.setState(q,dq);
   myANEA = myANEA.setD(d);
   myANEA = myANEA.setDprior(d);
   myANEA = myANEA.setXprior([q; dq]);
   myANEA = myANEA.setXvariance(Sx);
   if norm(myANEA.dDb_s.matrix - dD) > 1e-3
      disp(['[DANEA] dD numerical derivative is quite different: ' num2str(norm(myDNEA.dDb_s.matrix - dD))])
      res = 1;
   end

   if norm(myANEA.dDb.matrix - dD) > 1e-3
      disp(['[DANEA] dD numerical derivative is quite different: ' num2str(norm(myDNEA.dDb.matrix - dD))])
      res = 1;
   end
   
   if norm(full(myANEA.Sx_inv.matrix) - inv(Sx)) ~= 0
      disp('[DANEA] Sx_inv was not properly set!')
      res = 1;
   end
   
   if norm(full(myANEA.Sx.matrix) - Sx) ~= 0
      disp('[DANEA] Sx was not properly set!')
      res = 1;
   end

end


