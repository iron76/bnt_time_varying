
%% CHECK RNEA 
% provides the inverse dynamics analysis assuming no external forces.
% It can be useful exploiting this analysis as a benchmark.
%
%% Toolbox requirements:
% - iDynTree - mex
% - Featherstone toolbox (v2) with ID corrected as in bnt_time_varying repository

clear;clc;close all;

%testOptions
plotJointAngles = false;
plotJointVelocities = false;
plotJointAccelerations = false;
analyseRNEA = false;


%% Load Drake model

load('humanThreeLinkModelFromURDF_subject1.mat');

%% Check model imported from Drake
% 
% leg_R_foot = humanThreeLink_dmodel.Xtree{1}(1:3,1:3);
% torso_R_leg = humanThreeLink_dmodel.Xtree{2}(1:3,1:3);
% leg_r_foot = skew(-leg_R_foot' * humanThreeLink_dmodel.Xtree{1}(4:6,1:3));


%% Load acquired data

load('preProcessedSensorData.mat','processedSensorData');

t = processedSensorData(1,1).t;
f_temp = processedSensorData(1,1).f_0;
[~,chosenF_ID] = max(f_temp(:,1));
[~,tminIndex] = max(f_temp(chosenF_ID,1:round(end/2)));
[val,tmaxIndex] = max(f_temp(chosenF_ID,round(end/2):end));
totPointsInConsideration = tmaxIndex - tminIndex;
tminIndex = tminIndex + round(0.1*totPointsInConsideration);
tmaxIndex = tmaxIndex + round(length(t)/2) - round(0.1*totPointsInConsideration);
data.min_time = t(tminIndex);
data.max_time = t(tmaxIndex);
data.nsamples = tmaxIndex - tminIndex;
data.time = t(tminIndex:tmaxIndex);


%% Getting variables

endIdx = 18607;
%GET JOINT ANGLES
q1 = processedSensorData.q1;
q1 = q1(tminIndex:tmaxIndex,:); %window filter
q2 = processedSensorData.q2;
q2 = q2(tminIndex:tmaxIndex,:); %window filter

%GET JOINT VELOCITIES 
dq1 = processedSensorData.dq1;
dq1 = dq1(tminIndex:tmaxIndex,:); %window filter
dq2 = processedSensorData.dq2;
dq2 = dq2(tminIndex:tmaxIndex,:); %window filter
% dq1 = zeros (length(q1),1);
% dq2 = zeros (length(q2),1);


%GET JOINT ACCELERATIONS
% ddq1 = processedSensorData.ddq1;
% ddq1 = ddq1(tminIndex:tmaxIndex,:); %window filter
% ddq2 = processedSensorData.ddq2;
% ddq2 = ddq2(tminIndex:tmaxIndex,:); %window filter
ddq1 = zeros (length(q1),1);
ddq2 = zeros (length(q2),1);

q  = [q1,q2];
dq  = [dq1,dq2];
ddq  = [ddq1,ddq2];

%% Plot joint angles

% fig = figure();
% axes1 = axes('Parent',fig,'FontSize',16);
% box(axes1,'on');
% hold(axes1,'on');
% grid on;
% 
% plot1 = plot(data.time,(180/pi)*q1,'lineWidth',3.0, 'Parent',axes1);hold on;
% set(plot1,'color',[1 0 0]);
% plot2 = plot(data.time,(180/pi)*q2,'lineWidth',3.0, 'Parent',axes1);hold on;
% set(plot2,'color',[0 0.498039215803146 0]);
% leg = legend('$q_1$','$q_2$','Location','northwest');
% set(leg,'Interpreter','latex');
% set(leg,'FontSize',18);
% xlabel('Time [s]','FontSize',20);
% ylabel('Joint Angle [deg]','FontSize',20);
% axis tight; 

%% Plot joint velocities

% fig = figure();
% axes1 = axes('Parent',fig,'FontSize',16);
% box(axes1,'on');
% hold(axes1,'on');
% grid on;
% 
% plot1 = plot(data.time,(180/pi)*dq1,'lineWidth',3.0, 'Parent',axes1);hold on;
% set(plot1,'color',[1 0 0]);
% plot2 = plot(data.time,(180/pi)*dq2,'lineWidth',3.0, 'Parent',axes1);hold on;
% set(plot2,'color',[0 0.498039215803146 0]);
% leg = legend('$\dot q_1$','$\dot q_2$','Location','northwest');
% set(leg,'Interpreter','latex');
% set(leg,'FontSize',18);
% xlabel('Time [s]','FontSize',20);
% ylabel('Joint Velocity [deg/sec]','FontSize',20);
% axis tight; 

%% Plot joint accelerations

% fig = figure();
% axes1 = axes('Parent',fig,'FontSize',16);
% box(axes1,'on');
% hold(axes1,'on');
% grid on;
% 
% plot1 = plot(data.time,(180/pi)*ddq1,'lineWidth',3.0, 'Parent',axes1);hold on;
% set(plot1,'color',[1 0 0]);
% plot2 = plot(data.time,(180/pi)*ddq2,'lineWidth',3.0, 'Parent',axes1);hold on;
% set(plot2,'color',[0 0.498039215803146 0]);
% leg = legend('$\ddot q_1$','$\ddot q_2$','Location','northwest');
% set(leg,'Interpreter','latex');
% set(leg,'FontSize',18);
% xlabel('Time [s]','FontSize',20);
% ylabel('Joint Acceleration [deg/s^2]','FontSize',20);
% axis tight; 

%% Computing tau using Newton-Euler with Featherstone toolbox

tau = zeros(size (q));
% a = cell (size(q));
% fB_i = cell (size(q));
% f_i = cell (size(q));

fB_1 = zeros(size(q,1),6);
a_2 = zeros(size(q,1),6);
for i = 1:size(q)
      [tau_i, a_i, fB_i, f_i] = ID( humanThreeLink_dmodel, q(i,:), dq(i,:), ddq(i,:));
      tau(i,:) = tau_i';
       a_2(i,:) = a_i{2}';
       fB_1(i,:) = fB_i{1}';
%       f(i,:) = f_i';   
end

%% raw ID components
figure(1);
subplot(311);
plot1 = plot(data.time,(180/pi)*q1,'lineWidth',2.0); hold on;
set(plot1,'color',[1 0 0]);
plot2= plot(data.time,(180/pi)*q2,'lineWidth',2.0); hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$q_1$','$q_2$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Joint Angle [deg]','FontSize',18);
axis tight;
grid on;
 
subplot(312);
plot1 = plot(data.time,(180/pi)*dq1,'lineWidth',2.0); hold on;
set(plot1,'color',[1 0 0]);
plot2= plot(data.time,(180/pi)*dq2,'lineWidth',2.0); hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\dot q_1$','$\dot q_2$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Joint Velocity [deg/s]','FontSize',18);
axis tight;
grid on;

subplot(313);
plot1 = plot(data.time,tau(:,1), 'lineWidth',2.0); hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,tau(:,2), 'lineWidth',2.0); hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\tau_1$','$\tau_2$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Torque [Nm]','FontSize',18);
axis tight;
grid on;


figure();
plot(data.time,a_2(:,4:6));
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('linear Spatial Acceleration of link 2 (m/sec^2)','FontSize',18);
legend('$a_x$','$a_y$','$a_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on;

figure();
subplot(211);
plot(data.time,fB_1(:,4:6), 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('ForcePlate Force (N)','FontSize',18);
title('Wrench transmitted from link 1 to base','FontSize',15);
        legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on;

subplot(212);
plot(data.time,fB_1(:,1:3), 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('ForcePlate Momment (Nm)','FontSize',18);
legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on;

%% Computing the sensor measurement in its frame

if(exist('transformsData.mat','file')~=0)
    load('transformsData.mat');
    X_imu_2 = transformsData.X_imu_2;
    XStar_fp_0 = transformsData.XStar_fp_0;
    XStar_fp_1 = transformsData.XStar_fp_1;
end
nT = size(a_2,1);

S_lin = [zeros(3) ones(3)];
a_imu = S_lin*(X_imu_2*[a_2']);

%f_fp = (XStar_0_fp *fB_1')';
f_fp = (XStar_fp_1 *fB_1')';

if(analyseRNEA == true)
    checkRNEA_iDynTree;
end
figure();
plot(data.time,a_imu, 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Accelerometer measurement (m/sec^2)','FontSize',18);
legend('$a_x$','$a_y$','$a_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on;

figure();
subplot(211);
plot(data.time,f_fp(:,4:6), 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('ForcePlate Force (N)','FontSize',18);
title('Wrench measured in forceplate','FontSize',15);
        legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on;

subplot(212);
plot(data.time,f_fp(:,1:3), 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('ForcePlate Momment (Nm)','FontSize',18);
legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on;

 %% Plotting results
save('./experiments/humanFixedBase/resultsFromCheckRNEA.mat', 'tau');