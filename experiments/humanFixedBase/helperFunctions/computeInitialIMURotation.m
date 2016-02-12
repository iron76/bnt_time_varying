function [ R_G_imu,P_G_imu ] = computeInitialIMURotation(P_G_imuAT,P_G_imuBT,P_G_imuCT)
%COMPUTEINITIALIMUROTATION Function computes the position rotation of the IMU frame in the
% Global frame G. The IMU has 3 markers, A, B, and C. The marker A coincides with the IMU origin.


% Standard configuration:
%
%   (B)
%    |                         ^ x
%    |                         |
%    |                         |
%   (C)-------(A)        y <---o  z (exits, origin in A)



    numValues = 10;
  	P_G_imuA = mean(P_G_imuAT(1:numValues,:),1);
    P_G_imuB = mean(P_G_imuBT(1:numValues,:),1);
    P_G_imuC = mean(P_G_imuCT(1:numValues,:),1);
   
    % assumption: (B-C) coincides with x axes of imu
    x_G_imu = P_G_imuB - P_G_imuC;
    x_G_imu = x_G_imu ./ norm(x_G_imu); %norm
    
    % since (A-C) may not coincide with y axes of imu (assumption) 
    y_G_imuCandidate = P_G_imuC - P_G_imuA ;
    y_G_imuCandidate = y_G_imuCandidate ./ norm(y_G_imuCandidate); %norm
    
    y_G_imu = y_G_imuCandidate - (dot(y_G_imuCandidate',x_G_imu')') .* x_G_imu;
    y_G_imu = y_G_imu ./ norm(y_G_imu); %norm
    
    z_G_imu = cross(x_G_imu,y_G_imu);
    
    R_G_imu = [x_G_imu',y_G_imu',z_G_imu'] ;
    
    P_G_imu = P_G_imuA; % frame of IMU placed in A

end

