function [ ymodel ] = iCubSens( dmodel )
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
% This specific function creates the model for the iCub sensors, which at
% the moment of writing the present function include: 
%
%  - F/T sensor in the right/left thigh
%  - F/T sensor in the right/left foot
%  - F/T sensor in the right/left bicep
%  - acc sensor in the right/left foot
%  - acc/gyro in the right/left hand
%  - acc/gyro in the head
%
% Author: Francesco Nori
% Genova, Dec 2014

ymodel.NB = dmodel.NB;
ny = 0;
            %6FT   %2ACC    %3ACC+GYRO
ymodel.ny =   6  +    2   +   3;
ymodel.m  = 6*6  +  3*2   + 6*3;
ymodel.Y  = cell(ymodel.ny, ymodel.NB);
ymodel.Ys = cell(ymodel.ny, ymodel.NB);


for i = 1 : dmodel.NB
   
   if strcmp(dmodel.linkname{i}, 'l_thigh') || strcmp(dmodel.linkname{i}, 'r_thigh')
      ny = ny + 1;
      dy = 6;
      ymodel.labels{ny,1} = strcat(dmodel.linkname{i}, '_FT');
      for j = 1 : dmodel.NB
         ymodel.Y{ny,j}   = zeros(dy,26);
         ymodel.Ys{ny,j}  = sparse(zeros(dy,26));
      end
      ymodel.Y{ny,i}      = [zeros(dy,6) zeros(dy,6) eye(dy) zeros(dy, 1) zeros(dy,6) zeros(dy, 1)];      
      ymodel.Ys{ny,i}     = sparse(1:6,13:18,ones(dy,1), dy, 26);
      ymodel.sizes{ny,1}  = dy;
   end
   
   if strcmp(dmodel.linkname{i}, 'l_upper_foot+l_foot+l_sole') || strcmp(dmodel.linkname{i}, 'r_upper_foot+r_foot+r_sole')
      ny = ny + 1;
      dy = 6;
      ymodel.labels{ny,1} = strcat(dmodel.linkname{i}, '_FT');
      for j = 1 : dmodel.NB
         ymodel.Y{ny,j}   = zeros(dy,26);         
         ymodel.Ys{ny,j}  = sparse(zeros(dy,26));
      end
      ymodel.Y{ny,i}      = [zeros(dy,6) zeros(dy,6) eye(dy) zeros(dy, 1) zeros(dy,6) zeros(dy, 1)];      
      ymodel.Ys{ny,i}     = sparse(1:6,13:18,ones(dy,1), dy, 26);
      ymodel.sizes{ny,1}  = dy;
   end

   if strcmp(dmodel.linkname{i}, 'l_upper_arm+l_arm') || strcmp(dmodel.linkname{i}, 'r_upper_arm+r_arm')
      ny = ny + 1;
      dy = 6;
      ymodel.labels{ny,1} = strcat(dmodel.linkname{i}, '_FT');
      for j = 1 : dmodel.NB
         ymodel.Y{ny,j}   = zeros(dy,26);         
         ymodel.Ys{ny,j}  = sparse(zeros(dy,26));
      end
      ymodel.Y{ny,i}      = [zeros(dy,6) zeros(dy,6) eye(dy) zeros(dy, 1) zeros(dy,6) zeros(dy, 1)];      
      ymodel.Ys{ny,i}     = sparse(1:6,13:18,ones(dy,1), dy, 26);
      ymodel.sizes{ny,1}  = dy;
   end
   
   if strcmp(dmodel.linkname{i}, 'l_upper_foot+l_foot+l_sole') || strcmp(dmodel.linkname{i}, 'r_upper_foot+r_foot+r_sole')
      ny = ny + 1;
      dy = 3;
      ymodel.labels{ny,1} = strcat(dmodel.linkname{i}, '_acc');
      for j = 1 : dmodel.NB
         ymodel.Y{ny,j}   = zeros(dy,26);         
         ymodel.Ys{ny,j}  = sparse(zeros(dy,26));
      end
      ymodel.Y{ny,i}      = [zeros(dy,3) eye(dy) zeros(dy,6) zeros(dy,6) zeros(dy, 1) zeros(dy,6) zeros(dy, 1)];      
      ymodel.Ys{ny,i}     = sparse(1:3,4:6,ones(dy,1), dy, 26);
      ymodel.sizes{ny,1}  = dy;
   end
   
   if strcmp(dmodel.linkname{i}, 'chest+torso+neck_1+neck_2+head+imu_frame')
      ny = ny + 1;
      dy = 6;
      ymodel.labels{ny,1} = strcat(dmodel.linkname{i}, '_acc+gyr');
      for j = 1 : dmodel.NB
         ymodel.Y{ny,j}   = zeros(dy,26);         
         ymodel.Ys{ny,j}  = sparse(zeros(dy,26));
      end
      ymodel.Y{ny,i}      = [eye(dy) zeros(dy,6) zeros(dy,6) zeros(dy, 1) zeros(dy,6) zeros(dy, 1)];      
      ymodel.Ys{ny,i}     = sparse(1:6,1:6,ones(dy,1), dy, 26);
      ymodel.sizes{ny,1}  = dy;
   end
      
   if strcmp(dmodel.linkname{i}, 'l_forearm+l_wrist_1+l_hand+l_gripper') || strcmp(dmodel.linkname{i}, 'r_forearm+r_wrist_1+r_hand+r_gripper')
      ny = ny + 1;
      dy = 6;
      ymodel.labels{ny,1} = strcat(dmodel.linkname{i}, '_acc+gyr');
      for j = 1 : dmodel.NB
         ymodel.Y{ny,j}   = zeros(dy,26);         
         ymodel.Ys{ny,j}  = sparse(zeros(dy,26));
      end
      ymodel.Y{ny,i}      = [eye(dy) zeros(dy,6) zeros(dy,6) zeros(dy, 1) zeros(dy,6) zeros(dy, 1)];      
      ymodel.Ys{ny,i}     = sparse(1:6,1:6,ones(dy,1), dy, 26);
      ymodel.sizes{ny,1}  = dy;
   end   
end

m  = sum(cell2mat(ymodel.sizes));

if ymodel.ny ~= ny || ymodel.m ~= m
   error('Something wrong with the number of sensors')
end

for i = 1 : dmodel.NB
   ymodel.ny = ymodel.ny + 1;
   ymodel.sizes{ymodel.ny,1} = 6;
   ymodel.labels{ymodel.ny,1} = ['fx'  num2str(i)];
   for j = 1 : dmodel.NB
      ymodel.Y{ymodel.ny,j}   = zeros(6,26);
      ymodel.Ys{ymodel.ny,j}  = sparse(zeros(6,26));
   end
   ymodel.Y{ymodel.ny,i}  = [zeros(6,6) zeros(6,6) zeros(6, 6) zeros(6,1) eye(6) zeros(6, 1)];
   ymodel.Ys{ymodel.ny,i} = sparse(1:6,20:25,ones(6,1), 6, 26);
end

ymodel.m  = sum(cell2mat(ymodel.sizes));

for i = 1 : ymodel.ny
  for j = 1 : ymodel.NB
    Yx{i,j}  = ymodel.Ys{i,j}(:, 1:19);
    Yy{i,j}  = ymodel.Ys{i,j}(:, 20:end);
  end
end

ymodel.Ys = cell2mat([Yx Yy]);
end