%% loading stored and interpolated data

load('IMU_VICON_ShiftedData.mat');
subjectID = 1;
trialID = 1 ; 

temp = imu_vicon_shiftedData(subjectID,trialID);


Fx = temp.f_fp(:,1);
Fy = temp.f_fp(:,2);
Fz = temp.f_fp(:,3);

Mx = temp.f_fp(:,1);
My = temp.f_fp(:,2);
Mz = temp.f_fp(:,3);

dz = 43.3; % mm for AMTI
P_FP_x = -(My+Fx*dz)./Fz;
P_FP_y =  (Mx-Fy*dz)./Fz;


figure(1);
plot(P_FP_x,P_FP_y);
hold on; 
plot(P_FP_x(1),P_FP_y(1),'ro'); axis tight;

figure(2);
plot(temp.t_vicon,P_FP_x,'r',temp.t_vicon,P_FP_y,'b'); axis tight;
legend('x','y');

%compute new COP_2, 
pSelec = size(temp.P_G_lhee,1);
P_G_lhee = temp.P_G_lhee(1:pSelec,:);
P_G_ltoe = temp.P_G_ltoe(1:pSelec,:);
P_G_rhee = temp.P_G_rhee(1:pSelec,:);
P_G_rtoe = temp.P_G_rtoe(1:pSelec,:);

[R_G_0,P_G_0] = computeFootRotation(P_G_lhee,P_G_rhee,P_G_ltoe,P_G_rtoe); 
disp(P_G_0)
