clear all
close all
clc

joint_to_test = {'Rx', 'Ry', 'Rz', 'Px', 'Py', 'Pz', 'R', 'P'};
% joint_to_test = {'Ry'};

num_of_tests  = 2;

for i = 1 : length(joint_to_test)
   for j = 1 : num_of_tests

      q = rand(1,1);
      f = @(x) jcalcderiv(joint_to_test{i}, x);
      
      dXjdq_num   = deriv(f, q);      
      [~,~,dXjdq] = jcalcderiv( joint_to_test{i}, q );
            
      if norm(dXjdq - dXjdq_num) > 1e-7
         disp(['dX numerical derivative is quite different: ' num2str(norm(dXjdq - dXjdq_num))])
      end
   end
end

run('iCub.m')
dmodel  = iCub_dmodel;
ymodel  = autoSensRNEA(dmodel);

dmodel  = autoTreeStochastic(dmodel, 1e-3);
ymodel  = autoSensStochastic(ymodel, 1e-2);

myModel = model(dmodel);
mySens  = sensors(ymodel);
mySNEA  = SNEA(myModel, mySens);

for i = 1 : num_of_tests
   for j = 1 : dmodel.NB
      q  = rand(dmodel.NB     ,1);
      dq = rand(dmodel.NB     ,1);
      d  = ones(dmodel.NB * 26,1);
      f  = @(x) computeV(mySNEA, x, j);
      dv = deriv(f, [q; dq]);
      
      mySNEA = mySNEA.setState(q,dq);
      if norm(cell2mat(mySNEA.dvdx(j,:)) - dv) > 1e-5
         disp(['dv numerical derivative is quite different: ' num2str(norm(cell2mat(mySNEA.dvdx(j,:)) - dv))])
      end
   end
end

for i = 1 : num_of_tests
   q  = rand(dmodel.NB     ,1);
   dq = rand(dmodel.NB     ,1);
   d  = ones(dmodel.NB * 26,1);
   f  = @(x) computeDb(mySNEA , d , x);
   dD = deriv(f, [q; dq]);
   
   mySNEA = mySNEA.setState(q,dq);
   mySNEA = mySNEA.setD(d);
   if norm(mySNEA.Ddbx - dD) > 1e-5
      disp(['dD numerical derivative is quite different: ' num2str(norm(mySNEA.Ddbx - dD))])
   end
   
end



