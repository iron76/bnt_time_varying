function [ ymodel ] = autoSensANEA( dmodel , use_a )
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
   ymodel.sizes{ymodel.ny,1} = 1;
   ymodel.labels{ymodel.ny,1} = ['d2q' num2str(i)];
   
   if (use_a)
      ymodel.ny = ymodel.ny + 1;
      ymodel.sizes{ymodel.ny,1} = 6;
      ymodel.labels{ymodel.ny,1} = ['a' num2str(i)];
   end
end

ymodel.m  = sum(cell2mat(ymodel.sizes));
ymodel.Y  = cell(ymodel.ny, ymodel.NB);
ymodel.Ys = cell(ymodel.ny, ymodel.NB);

for i = 1 : ymodel.ny
   for j = 1 : ymodel.NB
      my = ymodel.sizes{i,1};      
      if strcmp(ymodel.labels{i,1}, ['a' num2str(j)])
         ymodel.Y{i,j}  = [eye(6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,1:6,ones(6,1), my, 7);
      elseif strcmp(ymodel.labels{i,1}, ['d2q' num2str(j)])
         ymodel.Y{i,j}  = [zeros(1,6) eye(1)];
         ymodel.Ys{i,j} = sparse(1,7,1,my,7);
      else
         ymodel.Ys{i,j} = sparse(my, 7);
         ymodel.Y{i,j}  = zeros(my, 7);
      end      
   end
end

for i = 1 : ymodel.ny
  for j = 1 : ymodel.NB
    Yx{i,j}  = ymodel.Ys{i,j}(:, 1:6);
    Yy{i,j}  = ymodel.Ys{i,j}(:, 7:end);
  end
end

ymodel.Ys = cell2mat([Yx Yy]);
end
