function [ ymodel ] = autoSensSIEfriction( dmodel )
%AUTOSENSRNEA Generates an RNEA sensor distribution for the given dynamic model.
%   This function generates a structure that contains all measurements
%   needed to perform RNEA inverse dynanmic compuations on the supplied
%   articulated rigid body (dmodel). The model is assumed to be in 
%   featherstone-like format. The struct has the following fields:
%
%       ny - the number of available sensors
%
%       NB - the number of rigid links in the articulated rigid body
%
%   labels - a (ny,1) cell array containing the type of sensor. Possible type of
%            sensors are listed in the following where the subscript i 
%            refers to the i-th link of the articulated rigid body.
%                     - a_i :  spatial accelration
%                     - fB_i:  net spatial force 
%                     - f_i:   net sptaial force with parent
%                     - tau_i: joint torque
%                     - fx_i:  external spatial force
%                     - d2q_i: joint acceleration
%
%    sizes - a (ny,1) cell array containing the size of each sensor
%
%        m - the lenght of measurement vetor y (sum of all sizes)
%
%        Y - a (ny,NB) cell array describing the link between the sensors 
%            available sensors and the vector of variables d
%            describing the dynamics of the articulated rigid body. In
%            particulat we have:
%                       Y{i,j} d_j = y_i
%            being y_i the i-th measurement and d_j the following vector of
%            dynamic variables associated to the articulated rigid body:
%                       d_j = [a_j, fB_j, f_j, tau_j, fx_j, d2q_j]
%
%       Ys - a sparse representation of the whole-matrix Y
%
%   This specifc function generates a measurement matrix which includes
%   only the sensors avaialble in the SIE robot. These sensors are the
%   joint accelerations and the 4-load cells at the feet. The load cells
%   give a 2D projection of the 3D wrench applied at the feet. In practice
%   the load cell measures fz, tau_x and tau_y. Therfore, one of the 
%   available load cells is represented as a 2D wrench (force and torque) 
%   exchanged at the root link (denoted f(1) in the code). The second load
%   cell is represented as a 2D external wrench exchanged at the last link 
%   of the chain (denoted fx(N_B) in the code). 
%
% Author: Francesco Nori
% Genova, Dec 2014

ymodel.NB = dmodel.NB;

if ~checkModel(dmodel)
   error('You should provide a featherstone-like mdoel')
end

if dmodel.NB < 4
    error('The funciton autoSensSIE requires NB >= 6')
end

ymodel.ny = 0;
for i = 1 : ymodel.NB-1
   ymodel.ny = ymodel.ny + 1;
   ymodel.sizes{ymodel.ny,1} = 6;
   ymodel.labels{ymodel.ny,1} = ['fx'  num2str(i)];
end

for i = 3 : ymodel.NB-2
   ymodel.ny = ymodel.ny + 1;
   ymodel.sizes{ymodel.ny,1} = 1;
   ymodel.labels{ymodel.ny,1} = ['d2q' num2str(i)];
end

ymodel.ny = ymodel.ny + 1;
ymodel.sizes{ymodel.ny,1} = 5;
ymodel.labels{ymodel.ny,1} = ['f' num2str(1) '_2D_friction'];

ymodel.ny = ymodel.ny + 1;
ymodel.sizes{ymodel.ny,1} = 5;
ymodel.labels{ymodel.ny,1} = ['fx' num2str(ymodel.NB) '_2D_friction'];

ymodel.m  = sum(cell2mat(ymodel.sizes));
ymodel.Y  = cell(ymodel.ny, ymodel.NB);
ymodel.Ys = cell(ymodel.ny, ymodel.NB);

for i = 1 : ymodel.ny
   for j = 1 : ymodel.NB
      my = ymodel.sizes{i,1};
      ymodel.Y{i,j}  = zeros(my, 26);
      ymodel.Ys{i,j} = zeros(my, 26);
      if strcmp(ymodel.labels{i,1}, ['a' num2str(j)])
         ymodel.Y{i,j}  = [eye(6) zeros(6,6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,1:6,ones(6,1), 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['fB' num2str(j)])
         ymodel.Y{i,j}  = [zeros(6,6) eye(6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,7:12,ones(6,1), 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['f' num2str(j)])
         ymodel.Y{i,j}  = [zeros(6,6) zeros(6,6) eye(6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,13:18,ones(6,1), 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['tau' num2str(j)])
         ymodel.Y{i,j}  = [zeros(1,6) zeros(1,6) zeros(1,6) eye(1, 1) zeros(1,6) zeros(1, 1)];
         ymodel.Ys{i,j} = sparse(1,19,1,1,26);
      end
      if strcmp(ymodel.labels{i,1}, ['fx' num2str(j)])
         ymodel.Y{i,j}  = [zeros(6,6) zeros(6,6) zeros(6, 6) zeros(6,1) eye(6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,20:25,ones(6,1), 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['d2q' num2str(j)])
         ymodel.Y{i,j}  = [zeros(1,6) zeros(1,6) zeros(1,6) zeros(1, 1) zeros(1,6) eye(1, 1)];
         ymodel.Ys{i,j} = sparse(1,26,1,1,26);
      end
      if strcmp(ymodel.labels{i,1}, ['f' num2str(j) '_2D'])
         S = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 1];
         ymodel.Y{i,j}  = [zeros(3,6) zeros(3,6)          S  zeros(3, 1) zeros(3,6) zeros(3, 1)];
         ymodel.Ys{i,j} = sparse(1:3,[13 14 18], ones(3,1), 3, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['fx' num2str(j) '_2D'])
         S = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 1];
         ymodel.Y{i,j}  = [zeros(3,6) zeros(3,6) zeros(3, 6) zeros(3,1)          S  zeros(3, 1)];
         ymodel.Ys{i,j} = sparse(1:3,[20 21 25], ones(3,1), 3, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['f' num2str(j) '_2D_friction'])
         S = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
         ymodel.Y{i,j}  = [zeros(5,6) zeros(5,6)          S  zeros(5, 1) zeros(5,6) zeros(5, 1)];
         ymodel.Ys{i,j} = sparse(1:5,[13 14 16:18], ones(5,1), 5, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['fx' num2str(j) '_2D_friction'])
         S = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
         ymodel.Y{i,j}  = [zeros(5,6) zeros(5,6) zeros(5, 6) zeros(5,1)          S  zeros(5, 1)];
         ymodel.Ys{i,j} = sparse(1:5,[20 21 23:25], ones(5,1), 5, 26);
      end                  
   end
end


for i = 1 : ymodel.ny
  for j = 1 : ymodel.NB
    Yx{i,j}  = ymodel.Ys{i,j}(:, 1:19);
    Yy{i,j}  = ymodel.Ys{i,j}(:, 20:end);
  end
end

ymodel.Ys = cell2mat([Yx Yy]);
end
