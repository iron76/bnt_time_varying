clear all
close all
clc

joint_to_test = {'Rx', 'Ry', 'Rz', 'Px', 'Py', 'Pz', 'R', 'P'};
% joint_to_test = {'Ry'};

num_of_tests  = 1000;

for i = 1 : length(joint_to_test)
   for j = 1 : num_of_tests

      q = rand(1,1);
      f = @(x) jcalcderiv(joint_to_test{i}, x);
      
      dXjdq_num   = deriv(f, q);      
      [~,~,dXjdq] = jcalcderiv( joint_to_test{i}, q );
            
      if norm(dXjdq - dXjdq_num) > 1e-7
         disp(['Numerical derivative is quite different: ' num2str(norm(dXjdq - dXjdq_num))])
      end
   end
end