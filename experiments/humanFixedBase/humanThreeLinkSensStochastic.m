
function [ model ] = humanThreeLinkSensStochastic( model )
%HUMANFIVELINKSENSSTOCHASTIC Summary of this function goes here
%   Detailed explanation goes here
%AUTOSENSSTOCHASTIC Add stochastic component to a sensor distribution
%   This function takes a structure containing sensor distribution (similar
%   to the one created by autoSensSNEA) and adds to the structure some
%   fields that are used to represent the variance of the measurements. The
%   sensor model is the following:
%
%                 Y(q,dq) d = y
%
%   and the variance is associated with the confidence on the measurement
%   equation.


idSy_inv = []; jdSy_inv = []; dSy_inv=[];
my = 1;

%uncertainties
% so       = 10e+4 * 0.0013;%10*0.7805; % IMU Gyroscope
% sa       = 10e+4 * 0.00426;%10*0.05; % IMU Acclerometer
% sf       = 10e+4 * 52.9097;%10; %FT - force
% su       = 10e+4 * 0.031; % FT - momment
% sd       = 10e+4 * 50*0.2; % joint acceleration

so       = 0.7805; %IMU accelerometer
sa       = 0.05; %IMU gyroscope
sf       = 10; % FT moment
su       = 0.1;  %0.031; % FT force
sd       = 0.2; % joint acceleration

imuS = [sa*eye(3) zeros(3,3); zeros(3,3) so*eye(3)];
ftsS = [sf*eye(3) zeros(3,3); zeros(3,3) su*eye(3)];
accS = imuS(1:3, 1:3);
ftxS = ftsS.*1e-4;
d2qS = sd;

for i = 1 : model.ny
   dy = model.sizes{i,1};
   if strcmp(model.labels{i}(end-2:end), 'imu')
      model.Sy{i,1} = imuS;
   elseif strcmp(model.labels{i}(end-2:end), 'fts')
      model.Sy{i,1} = ftsS;
   elseif strcmp(model.labels{i}(end-2:end), 'acc')
      model.Sy{i,1} = accS;
   elseif strcmp(model.labels{i}(end-2:end), 'ftx')
      model.Sy{i,1} = ftxS;
   elseif strcmp(model.labels{i}(end-2:end), 'd2q')
      model.Sy{i,1} = d2qS;      
   end
      
      
   [ii, jj, ss] = placeSubmatrixSparse(my, my, inv(model.Sy{i,1}));
   idSy_inv = [idSy_inv; ii];
   jdSy_inv = [jdSy_inv; jj];
   dSy_inv  = [dSy_inv;  ss];
   my = my + dy;
end
model.Sy_inv = sparse(idSy_inv, jdSy_inv, dSy_inv);
end