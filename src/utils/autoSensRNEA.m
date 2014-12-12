function [ ymodel ] = autoSensRNEA( dmodel )
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
% Author: Francesco Nori
% Genova, Dec 2014

ymodel.NB = dmodel.NB;

if ~checkModel(dmodel)
   error('You should provide a featherstone-like mdoel')
end

ymodel.ny = 0;
for i = 1 : ymodel.NB
   ymodel.ny = ymodel.ny + 1;
   ymodel.sizes{ymodel.ny,1} = 6;
   ymodel.labels{ymodel.ny,1} = ['fx'  num2str(i)];
   ymodel.ny = ymodel.ny + 1;
   ymodel.sizes{ymodel.ny,1} = 1;
   ymodel.labels{ymodel.ny,1} = ['d2q' num2str(i)];
end

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
