subjectID = 3;
trialID = 1;

%% load dataset
load(sprintf('./experiments/humanFixedBase/savedBERDYresult_subj%d_trial%d.mat',subjectID,trialID),'res','data','myPNEA');

%% plot joint angles
figure(1);
plot(data.time,(180/pi)*data.q1,'r',data.time,(180/pi)*data.q2,'b','lineWidth',2.0);
legend('\theta_1','\theta_2');
xlabel('Time t(sec)');
ylabel('Joint Angle (degs)');

axis tight; 

%% plot joint velocities
figure(2);
plot(data.time,(180/pi)*data.dq1,'r',data.time,(180/pi)*data.dq2,'b','lineWidth',2.0);
legend('\theta_1','\theta_2');
xlabel('Time t(sec)');
ylabel('Joint Velocity (degs/sec)');
axis tight; 

%% plot IMU data
[RNEA_imu,RNEA_ftx,RNEA_t] = sensorDataConsistencyCheck(subjectID,trialID,'noplots');

figure(3);
subplot(2,1,1);
plot(data.time,data.y_imu(1:3,:)','lineWidth',1.0);hold on;
plot(RNEA_t,RNEA_imu(:,1:3),'--');
legend('a_x','a_y','a_z');
xlabel('Time t(sec)');
ylabel('Acceleration (m/sec^2)');
axis tight; 

subplot(2,1,2);
plot(data.time,data.y_imu(4:6,:)','lineWidth',1.0);hold on;
plot(RNEA_t,RNEA_imu(:,4:6),'--');
legend('\omega_x','\omega_y','\omega_z');
xlabel('Time t(sec)');
ylabel('Angular Velocity (rad/sec)');
axis tight; 

%% plot forcePlate

figure(4);
subplot(2,1,1);
plot(data.time,data.y_ftx(1:3,:)','lineWidth',1.0);hold on;
plot(RNEA_t,RNEA_ftx(:,1:3),'--');
legend('f_x','f_y','f_z');
xlabel('Time t(sec)');
ylabel('Force (N)');
axis tight; 

subplot(2,1,2);
plot(data.time,data.y_ftx(4:6,:)','lineWidth',1.0);hold on;
plot(RNEA_t,RNEA_ftx(:,4:6),'--');
legend('\mu_x','\mu_y','\mu_z');
xlabel('Time t(sec)');
ylabel('Momment (N/m)');
axis tight; 


%% plot torque estimates

figure(5);
plot(data.time,res.tau_foot,'b',data.time,res.tau_leg,'r','lineWidth',2.0);
xlabel('Time (sec)');
ylabel('Torques (Nm)');
axis tight; hold on;   
legend('\tau_1 (ankle)','\tau_2 (hip)');


%% plot torque variance
figure(6);
shadedErrorBar(data.time,res.tau_foot,squeeze(2*sqrt(res.Stau_foot)),'b');hold on;
shadedErrorBar(data.time,res.tau_leg,squeeze(2*sqrt(res.Stau_foot)),'r');

%plot(tempT,res.tau_foot,'b',tempT,res.tau_leg,'r','lineWidth',2.0);
xlabel('Time (sec)');
ylabel('Torques (Nm)');
axis tight; hold on;   
legend('\tau_1 (ankle)','\tau_2 (hip)');


%% save plots
figBaseFolder = './experiments/humanFixedBase/plots';
if(exist(figBaseFolder,'dir')==0)
    mkdir(figBaseFolder);
end

plotFigBaseName = strcat(figBaseFolder,sprintf('/PNEA_human%d_trial%d_',subjectID,trialID));
selectedFigList = [1,2,3,4,5,6];
FigName = {'jointAngles',...
            'jointVelocities',...
            'imu',...
            'forcePlate',...
            'torques',...
            'TorquesVariance'};

        
for i = 1:length(selectedFigList)
     figure(selectedFigList(i))
     set(gca,'FontSize',12);
     set(gcf,'Renderer','OpenGL');
     print('-depsc2','-r200',strcat(plotFigBaseName,FigName{i}),'-opengl');
end