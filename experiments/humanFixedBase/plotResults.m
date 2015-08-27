load('./experiments/humanFixedBase/savedBERDYresult.mat');


tempT = data.time;%data.time(1:10:end-3);
figure;
plot(tempT,res.tau_foot,'b','LineWidth',2); hold on
plot(tempT,res.tau_leg,'r','LineWidth',2);
xlabel('Time (sec)');
ylabel('Torques (Nm)');
axis tight; hold on;   
legend('\tau_1 (ankle)','\tau_2 (hip)','Location','NorthEast');
title('Torque Estimates');


tempT = data.time;%data.time(1:10:end-3);
figure;plot(tempT,data.q1,'b','LineWidth',2); hold on;
plot(tempT,data.q2,'r','LineWidth',2);
xlabel('Time (sec)');
ylabel('Joint Position (Deg)');
axis tight; hold on;   
legend('\theta_1 (ankle)','\theta_2 (hip)');
title('Joint Angles measured');

figure;

subplot(2,1,1);plot(tempT,data.y_imu(1:3,:),'LineWidth',2);
xlabel('Time (sec)');
ylabel('Acceleration (m/sec^2)');
axis tight; hold on;
title('IMU sensor readings');
legend('a_x','a_y','a_z');
subplot(2,1,2);plot(tempT,data.y_imu(4:6,:),'LineWidth',2);
xlabel('Time (sec)');
ylabel('Angular Velocity  (rad/sec)');
axis tight; hold on;
legend('\omega_x','\omega_y','\omega_z');
%legend('\theta_1 (ankle)','\theta_2 (hip)');


figure;
subplot(2,1,1);
plot(tempT,data.y_ftx(1:3,:),'LineWidth',2);
xlabel('Time (sec)');
ylabel('Force (N)');
axis tight; hold on;   
legend('F_x','F_y','F_z');
%legend('\theta_1 (ankle)','\theta_2 (hip)');
title('Force Plate sensor readings');

subplot(2,1,2);
plot(tempT,data.y_ftx(4:6,:),'LineWidth',2);
xlabel('Time (sec)');
ylabel('Torques (Nm)');
legend('\mu_x','\mu_y','\mu_z');
axis tight; hold on;   

