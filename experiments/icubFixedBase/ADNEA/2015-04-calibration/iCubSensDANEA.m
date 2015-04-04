function [ ymodel ] = iCubSensDNEA( dmodel , ymodel_ini, sens, mask_q, mask_dq )
%AUTOSENSSNEA Generates a random sensor distribution articulated rigid body.
%   This function generates a structure that contains all measurements
%   needed to perform sparse inverse dynanmic compuations on the supplied
%   articulated rigid body (dmodel) assuming mutiple (and possibly
%   redudant) senors. The dynamic model is assumed to be in
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
%                     - dq_i:  joint velocity
%                     -  q_i:  joint position
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
%                       d_j = [a_j, fB_j, f_j, tau_j, fx_j, d2q_j; q_j; dq_j]
%
%       Ys - a sparse representation of the whole-matrix Y
%
% Author: Francesco Nori
% Genova, Dec 2014

ymodel = ymodel_ini;

ymodel.Ys  = [ymodel.Ys, sparse([],[],[],ymodel.m, 2*dmodel.NB,0)];
for i = 1 : ymodel.ny
   ymodel.Y{i,dmodel.NB+1} = zeros(size(ymodel.Y{i,1},1), dmodel.NB);
   ymodel.Y{i,dmodel.NB+2} = zeros(size(ymodel.Y{i,1},1), dmodel.NB);
end

%% Angular velocity measurements
for i = 1 : length(sens.parts)
   sens_label = sens.labels{i};
   if strcmp(sens_label(end-2:end), 'gyr')
      
      for j = 1 : ymodel.NB
         if(strcmp(sens.parts{i}, dmodel.linkname{j}))
            link_ind = j;
         end
      end
      
      ymodel.ny = ymodel.ny + 1;
      ymodel.sizes{ymodel.ny,1} = 3;
      ymodel.m  = ymodel.m  + ymodel.sizes{ymodel.ny,1};
      ymodel.labels{ymodel.ny,1} = ['y_omega' num2str(link_ind)];
      
      for j = 1 : dmodel.NB
         ymodel.Y{ymodel.ny,j} = zeros(ymodel.sizes{ymodel.ny,1}, 7);
      end
      ymodel.Y{ymodel.ny,dmodel.NB+1} = zeros(ymodel.sizes{ymodel.ny,1}, dmodel.NB);
      ymodel.Y{ymodel.ny,dmodel.NB+2} = zeros(ymodel.sizes{ymodel.ny,1}, dmodel.NB);
      ymodel.Ys = [ ymodel.Ys;  sparse(zeros(ymodel.sizes{ymodel.ny,1}, 9*dmodel.NB)) ];
   end
end

%% [q; dq]
% for i = 1 : length(mask_q)
%    if mask_q(i) == 1
%       yi = zeros(1, 2*dmodel.NB);
%       yi(1,i) = 1;
%       ymodel.Ys =      [ ymodel.Ys;  [sparse(zeros(1, 26*dmodel.NB)) sparse(yi)]];
%       for j = 1 : dmodel.NB
%          ymodel.Y{ymodel.ny+1,j} = zeros(1, 26);
%       end
%       ymodel.Y{ymodel.ny+1,dmodel.NB+1} = yi(1, 1:dmodel.NB);
%       ymodel.Y{ymodel.ny+1,dmodel.NB+2} = yi(1, dmodel.NB+1:end);
%       ymodel.m  = ymodel.m  + 1;
%       ymodel.ny = ymodel.ny + 1;
%       ymodel.sizes{ymodel.ny,1} = 1;
%       ymodel.labels{ymodel.ny,1} = ['y_q' num2str(i)];
%    end
% end
for i = 1 : length(mask_dq)
   if mask_dq(i) == 1
      yi = zeros(1, 2*dmodel.NB);
      yi(1, dmodel.NB + i) = 1;
      ymodel.Ys = [ ymodel.Ys;  sparse(zeros(1, 7*dmodel.NB)) sparse(yi)];
      for j = 1 : dmodel.NB
         ymodel.Y{ymodel.ny+1,j} = zeros(1, 7);
      end
      ymodel.Y{ymodel.ny+1,dmodel.NB+1} = yi(1, 1:dmodel.NB);
      ymodel.Y{ymodel.ny+1,dmodel.NB+2} = yi(1, dmodel.NB+1:end);
      ymodel.m  = ymodel.m  + 1;
      ymodel.ny = ymodel.ny + 1;
      ymodel.sizes{ymodel.ny,1} = 1;
      ymodel.labels{ymodel.ny,1} = ['y_dq' num2str(i)];
   end
end
