function res = testDerivatives(dmodel, ymodel)

res = 0;

joint_to_test = {'Rx', 'Ry', 'Rz', 'Px', 'Py', 'Pz', 'R', 'P'};
% joint_to_test = {'Ry'};

num_of_tests  = 1;

for i = 1 : length(joint_to_test)
   for j = 1 : num_of_tests

      q = rand(1,1);
      f = @(x) jcalcderiv(joint_to_test{i}, x);
      
      dXjdq_num   = deriv(f, q);      
      [~,~,dXjdq] = jcalcderiv( joint_to_test{i}, q );
            
      if norm(dXjdq - dXjdq_num) > 1e-7
         disp(['[DERIVATIVES] dX numerical derivative is quite different: ' num2str(norm(dXjdq - dXjdq_num))])
         res = 1;
      end
   end
end

myModel = model(dmodel);
mySens  = sensors(ymodel);
myDNEA  = DNEA(myModel, mySens);

% for i = 1 : num_of_tests
%    q  = rand(dmodel.NB     ,1);
%    dq = rand(dmodel.NB     ,1);
% 
%    myDNEA = myDNEA.setState(q,dq);   
%    for j = 1 : NB
%       f  = @(x) computeV(myDNEA, x, j);
%       dv = deriv(f, [q; dq]);
%       
%       if norm(cell2mat(myDNEA.dvdx(j,:)) - dv) > 1e-5
%          disp(['dv numerical derivative is quite different: ' num2str(norm(cell2mat(myDNEA.dvdx(j,:)) - dv))])
%          res = 1;
%       end
%    end
% end

for i = 1 : num_of_tests
   q  = rand(dmodel.NB     ,1);
   dq = rand(dmodel.NB     ,1);
   d  = rand(dmodel.NB * 26,1);
   
   myDNEA = myDNEA.setState(q,dq);
   
   f  = @(x) computeY(myDNEA, d, x);
   dy = deriv(f, [q; dq]);
   
   % In current formulation 
   %
   % y = Yd d + Yx x + b_Y(x) = [Yd(x) Yx][d ; x] + b_Y(x)
   %
   % and therefore the derivative is:
   %
   % dy/dx = Yx + db_Y/dx
   
   Yx  = myDNEA.IDsens.sensorsParams.Ys(:,(end-2*dmodel.NB+1):end);
   dby = myDNEA.dby_s.matrix + Yx;
   if norm( dby - dy) > 1e-6
      disp(['[DERIVATIVES] dy numerical derivative is quite different: ' num2str(norm(dby - dy))])
      imagesc([dby - dy])
      colorbar
      res = 1;
   end
end


for i = 1 : num_of_tests
   q  = rand(dmodel.NB     ,1);
   dq = rand(dmodel.NB     ,1);
   d  = rand(dmodel.NB * 26,1);
   Sx = diag(rand(dmodel.NB*2, 1));
   f  = @(x) computeDb(myDNEA , d , x);
   dD = deriv(f, [q; dq]);
   
   myDNEA = myDNEA.setState(q,dq);
   myDNEA = myDNEA.setD(d);
   myDNEA = myDNEA.setDprior(d);
   myDNEA = myDNEA.setXprior([q; dq]);
   myDNEA = myDNEA.setXvariance(Sx);
   myDNEA.solveID();
   if norm(myDNEA.dDb_s.matrix - dD) > 1e-3
      disp(['[DERIVATIVES] dD numerical derivative is quite different: ' num2str(norm(myDNEA.dDb_s.matrix - dD))])
      res = 1;
   end

   if norm(myDNEA.dDb.matrix - dD) > 1e-3
      disp(['[DERIVATIVES] dD numerical derivative is quite different: ' num2str(norm(myDNEA.dDb.matrix - dD))])
      res = 1;
   end
   
   if norm(full(myDNEA.Sx_inv.matrix) - inv(Sx)) ~= 0
      disp('[DERIVATIVES] Sx_inv was not properly set!')
      res = 1;
   end
   
   if norm(full(myDNEA.Sx.matrix) - Sx) ~= 0
      disp('[DERIVATIVES] Sx was not properly set!')
      res = 1;
   end

end



