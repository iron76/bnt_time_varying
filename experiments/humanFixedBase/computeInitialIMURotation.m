function [ R_G_imu0,P_G_imu0 ] = computeInitialIMURotation(P_G_imuAT,P_G_imuBT,P_G_imuCT)
%COMPUTEINITIALIMUROTATION Summary of this function goes here
%   Detailed explanation goes here
    
    %considering initial values only
    numValues = 10;
  	P_G_imuA = mean(P_G_imuAT(1:numValues,:),1);
    P_G_imuB = mean(P_G_imuBT(1:numValues,:),1);
    P_G_imuC = mean(P_G_imuCT(1:numValues,:),1);
   
    x_G_imu = computeVectorFromPoints(P_G_imuC,P_G_imuA);
    y_G_imuPossible = computeVectorFromPoints(P_G_imuC,P_G_imuB);
    
    %normalising
    x_G_imu = x_G_imu ./ norm(x_G_imu);
    y_G_imuPossible = y_G_imuPossible./norm(y_G_imuPossible);
    y_G_imu = y_G_imuPossible - dot(x_G_imu',y_G_imuPossible')';
    
    z_G_imu = cross(x_G_imu,y_G_imu);
    R_G_imu0 = [x_G_imu',y_G_imu',z_G_imu'] ;
    P_G_imu0 = computeCentroidOfTriangle(P_G_imuA,P_G_imuB,P_G_imuC);
end

