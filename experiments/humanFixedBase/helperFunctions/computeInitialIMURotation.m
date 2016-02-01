function [ R_G_imu0,P_G_imu0 ] = computeInitialIMURotation(P_G_imuAT,P_G_imuBT,P_G_imuCT)
%COMPUTEINITIALIMUROTATION Function computes the position rotation of the IMU frame in the
% Global frame G. The IMU rotation is such that in the global frame, the
% IMU_x is along G_z, IMU_y is along -G_x and IMU_z is along -G_y.
% The IMU has 3 markers, A, B, and C. The B is higher than the other two
% when in the global frame. The marker A coincides with the IMU origin.
%
% Author: Naveen Kuppuswamy (naveen.kuppuswamy@iit.it)
% iCub Facility, Istituto Italiano di Tecnologia, 21 January 2016

    %considering initial values only
    numValues = 10;
  	P_G_imuA = mean(P_G_imuAT(1:numValues,:),1);
    
    
    P_G_imuB = mean(P_G_imuBT(1:numValues,:),1);
    P_G_imuC = mean(P_G_imuCT(1:numValues,:),1);
   
   P_G_imuFakeB = [P_G_imuA(1) P_G_imuA(2) P_G_imuB(3)];
    x_G_imu = computeVectorFromPoints(P_G_imuA,P_G_imuFakeB);
    x_G_imu_hat = x_G_imu ./ norm(x_G_imu);
    
    y_G_imuCandidate = computeVectorFromPoints(P_G_imuA,P_G_imuC);
    y_G_imuCandidate_hat = y_G_imuCandidate ./ norm(y_G_imuCandidate);
    
    y_G_imu = y_G_imuCandidate_hat - (dot(y_G_imuCandidate_hat',x_G_imu_hat')') .* x_G_imu_hat;
    y_G_imu_hat = y_G_imu ./ norm(y_G_imu);
    
    z_G_imu = cross(x_G_imu_hat,y_G_imu_hat);
    
    R_G_imu0 = [x_G_imu_hat',y_G_imu_hat',z_G_imu'] ;
    
    P_G_imu0 = P_G_imuA;
   
    
end

