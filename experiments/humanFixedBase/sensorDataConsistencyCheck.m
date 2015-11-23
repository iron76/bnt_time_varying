function [RNEA_imu,RNEA_ftx,t] = sensorDataConsistencyCheck(subjectID,trialID,plots) 
%subjectID = 3;
%trialID = 1;
data = struct();

%% load subject specific URDF
humanThreeLink
data = organiseBERDYCompatibleSensorData(data,subjectID,trialID);
%load('./experiments/humanFixedBase/organiseBerdyCompatibleSensorData.mat');
load('./experiments/humanFixedBase/data/preProcessedSensorData.mat','processedSensorData');
sensor_imu = data.y_imu';%processedSensorData(subjectID,trialID).imu;
sensor_ftx = data.y_ftx';%processedSensorData(subjectID,trialID).ftx;
R_G_2 = processedSensorData(subjectID,trialID).R_G_2;

q1 = data.q1;
q2 = data.q2;

dq1 = data.dq1;%zeros (length(q1),1);
dq2 = data.dq2;%zeros (length(q2),1);

ddq1 = data.ddq1;%zeros (length(q1),1);
ddq2 = data.ddq2;%zeros (length(q2),1);

t = data.time;

q  = [q1,q2];
dq  = [dq1,dq2];
ddq  = [ddq1,ddq2];


tau = zeros(size(q));
%a_1 = zeros(size(q,1),6);
a_2_2 = zeros(size(q,1),3);
%v_1 = zeros(size(q,1),6);
v_2_2 = zeros(size(q,1),3);
%fB_1 = zeros(size(q,1),6);
fBase = zeros(size(q,1),6);
%f_1 = zeros(size(q,1),6);
%f_2 = zeros(size(q,1),6);

R_D_G = eye(3);%[ 0 -1 0; 0 0 -1; 1 0 0];
%tau = zeros(size(q,1),2);

humanThreeLink_dmodel.gravity = [+0;9.81;0];
%a = cell (size(q));
for i = 1:size(q)
     ftx_sens{1} = sensor_ftx(i,:)';
     ftx_sens{2} = zeros(6,1);
     %q(i,:) = zeros(1,2);
     %dq(i,:) = zeros(1,2);
     %ddq(i,:) = zeros(1,2);
      [tau_i, a_i, fB_i, f_i, v_i,fBase_i] = ID( humanThreeLink_dmodel, q(i,:), dq(i,:), ddq(i,:));
      tau(i,:) = tau_i';      
      %a_1(i,:) = a_i{1}';
      a_2_2(i,:) = a_i{2}(4:6)';% + [0 0 -9.8];
      %v_1(i,:) = v_i{1}';
      v_2_2(i,:) = v_i{2}(1:3)';
      %fB_1(i,:)= fB_i{1}';
      fBase(i,:)= fBase_i';
      %f_1(i,:)= f_i{1}';
      %f_G_2(i,:)= (R_G_2{i}*f_i{2}(1:3))';
end

RNEA_imu = [a_2_2 v_2_2 ];
RNEA_ftx = [fBase(:,4:6) fBase(:,1:3)];

% 
% figure;
% plot(t,tau(:,1)); hold on;
% xlabel('time (sec)');
% ylabel('torque (Nm)');
% legend('\tau_1','\tau_2');
% axis tight;

tau2 = zeros(size(q));
for i = 1:size(q)
     ftx_sens{1} = sensor_ftx(i,:)';
     ftx_sens{2} = zeros(6,1);
      [tau2_i, a_i, fB_i, f_i, v_i,fBase_i] = ID( humanThreeLink_dmodel, q(i,:),0.* dq(i,:), 0.*ddq(i,:));
      tau2(i,:) = tau2_i';      
      %a_1(i,:) = a_i{1}';
      a_2_2(i,:) = a_i{2}(4:6)';% + [0 0 -9.8];
      %v_1(i,:) = v_i{1}';
      v_2_2(i,:) = v_i{2}(1:3)';
      %fB_1(i,:)= fB_i{1}';
      fBase(i,:)= fBase_i';
      %f_1(i,:)= f_i{1}';
      %f_G_2(i,:)= (R_G_2{i}*f_i{2}(1:3))';
end

if(strcmp(plots,'noplots')~=1)


    figure(1);
    subplot(2,1,1);
    plot(t,a_2_2);
    title('Accelerometer');
    xlabel('time (sec)');
    ylabel('predicted a_2 (m/sec^2)');
    legend('a_x','a_y','a_z');

    axis tight;

    subplot(2,1,2);
    plot(t,(R_D_G*sensor_imu(:,1:3)')');
    xlabel('time (sec)');
    ylabel('actual a_2 (m/sec^2)');
    legend('a_x','a_y','a_z');
    axis tight;


    figure(2);
    subplot(2,1,1);
    plot(t,v_2_2);
    title('Gyroscope');
    xlabel('time (sec)');
    ylabel('predicted \omega_2 (rad/sec)');
    legend('\omega_x','\omega_y','\omega_z');
    axis tight;

    subplot(2,1,2);
    plot(t,(R_D_G*sensor_imu(:,4:6)')');
    xlabel('time (sec)');
    ylabel('actual \omega_2 (rad/sec)');
    legend('\omega_x','\omega_y','\omega_z');
    axis tight;

    figure(3);
    subplot(2,1,1);
    plot(t,fBase(:,4:6));
    title('Base Force');
    xlabel('time (sec)');
    ylabel('predicted fBase (N)');
    legend('f_x','f_y','f_z');
    axis tight;

    subplot(2,1,2);
    plot(t,sensor_ftx(:,1:3));
    xlabel('time (sec)');
    ylabel('actual fBase (N)');
    legend('f_x','f_y','f_z');
    axis tight;
    
    figure(4);
    subplot(2,1,1);
    plot(t,fBase(:,1:3));
    title('Base Momment');
    xlabel('time (sec)');
    ylabel('predicted MuBase (Nm)');
    legend('\mu_x','\mu_y','\mu_z');
    axis tight;

    subplot(2,1,2);
    plot(t,sensor_ftx(:,4:6));
    xlabel('time (sec)');
    ylabel('actual fBase (N)');
    legend('\mu_x','\mu_y','\mu_z');
    axis tight;

    figure(5);
    plot(t,tau);
    xlabel('time (sec)');
    ylabel('Torrque (Nm)');
    legend('\tau_1','\tau_2');
    axis tight;
  
end
end