function [ ymodel ] = basicSensing( dmodel  )
%AUTOSENSSNEA Generates a sensor distribution articulated rigid body.
%   This function generates a structure that contains just: 
%           * Six Axis F/T sensor at the base 
%           * Joint accelerations 
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
%    If a second argument sens is specified in the form {'s1', 's2'} the
%    sensors 's1' and 's2' are included among sensors where possible
%    sensors are 'a', 'f', 'fB', 'd2q', 'tau', 'fB'.
%
% Author: Francesco Nori
% Genova, Dec 2014

if nargin == 1
   sens = [];
end

ymodel.NB = dmodel.NB;
ny = 0;

ny = ny + 1;
ymodel.sizes{ny,1} = 6;
ymodel.labels{ny,1} = ['y_f' num2str(1)];

for i = 1 : ymodel.NB
    ny = ny + 1;
    ymodel.sizes{ny,1} = 1;
    ymodel.labels{ny,1} = ['y_d2q' num2str(i)];
end

for i = 1 : ymodel.NB
    ny = ny + 1;
    ymodel.sizes{ny,1} = 6;
    ymodel.labels{ny,1} = ['y_fx' num2str(i)];
end

for i = 1 : ymodel.NB
    ny = ny + 1;
    ymodel.sizes{ny,1} = 6;
    ymodel.labels{ny,1} = ['y_a' num2str(i)];
end

ymodel.m  = sum(cell2mat(ymodel.sizes));
ymodel.ny = ny;
ymodel.Y  = cell(ny, ymodel.NB);
ymodel.Ys = cell(ny, ymodel.NB);

for i = 1 : ymodel.ny
   for j = 1 : ymodel.NB
      my = ymodel.sizes{i,1};
      ymodel.Y{i,j}  = zeros(my, 26);
      ymodel.Ys{i,j} = zeros(my, 26);
      if strcmp(ymodel.labels{i,1}, ['y_a' num2str(j)])
         d = ones(6,1);
         ymodel.Y{i,j}  = [diag(d) zeros(6,6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,1:6,d, 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['y_fB' num2str(j)])
         d = ones(6,1);
         ymodel.Y{i,j}  = [zeros(6,6) diag(d) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,7:12,d, 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['y_f' num2str(j)])
         d = ones(6,1);
         ymodel.Y{i,j}  = [zeros(6,6) zeros(6,6) diag(d) zeros(6, 1) zeros(6,6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,13:18,d, 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['y_tau' num2str(j)])
         d = ones(1,1);
         ymodel.Y{i,j}  = [zeros(1,6) zeros(1,6) zeros(1,6) diag(d) zeros(1,6) zeros(1, 1)];
         ymodel.Ys{i,j} = sparse(1,19,d,1,26);
      end
      if strcmp(ymodel.labels{i,1}, ['y_fx' num2str(j)])
         d = ones(6,1);
         ymodel.Y{i,j}  = [zeros(6,6) zeros(6,6) zeros(6, 6) zeros(6,1) diag(d) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,20:25,d, 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['y_d2q' num2str(j)])
         d = ones(1,1);
         ymodel.Y{i,j}  = [zeros(1,6) zeros(1,6) zeros(1,6) zeros(1, 1) zeros(1,6) diag(d)];
         ymodel.Ys{i,j} = sparse(1,26,d,1,26);
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
