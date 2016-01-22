% predictSensorMeasurement
% Script makes a prediction of the sensor measurements of the Human
% Dynamics Estimation Experiment. The sensors are the VICON markers, IMU
% placed on the chest and force place on the bottom of my foot. The method
% loads the sensor measurements (q, qDot) and the link sensor Transforms 
% from the intermediateDataFiles folder. 
%
% It can be useful exploiting this analysis as a benchmark.
%
% Toolbox requirements:
% - iDynTree - mex
% - Featherstone toolbox (v2) with ID corrected as in bnt_time_varying repository
% 
% Prerequisite : 
% Run synchroniseCaptureData, computeLinkSensorFrames and
% organiseSensorData if running this function for the first time
%
% Author: Naveen Kuppuswamy (naveen.kuppuswamy@iit.it)
% iCub Facility, Istituto Italiano di Tecnologia, 21 January 2016

clear;clc;close all;

%testOptions
plotJointAngles = false;
plotJointVelocities = false;
plotJointAccelerations = false;
analyseRNEA = false;


%% Load Drake model

load('humanThreeLinkModelFromURDF_subject1.mat');

%% Load acquired data

load('./experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat','processedSensorData');

t = processedSensorData(1,1).t;
f_measured = processedSensorData(1,1).f_fp;
a_measured = processedSensorData(1,1).a_imu;
tminIndex = 1;
tmaxIndex = length(t);
data.min_time = t(1);
data.max_time = t(end);
data.nsamples = length(t);
data.time = t;


%% Getting variables

endIdx = 1;%18607;
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

%% Computing tau using Newton-Euler with Featherstone toolbox

tau = zeros(size (q));
fB_1 = zeros(size(q,1),6);
a_2 = zeros(size(q,1),6);

fprintf('Iterating through time, computing dynamics using RNEA\n');

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

fprintf('Predicting sensor measurement using current sensor-link transforms\n');
%% Computing the sensor measurement in its frame

load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');
X_imu_2 = sensorLinkTransforms.X_imu_2;
XStar_fp_0 = sensorLinkTransforms.XStar_fp_0;
XStar_fp_1 = sensorLinkTransforms.XStar_fp_1;
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
title('Predicted vs Actual Accelerometer measurements');
legend('$a_x$','$a_y$','$a_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
hold on;
plot(data.time,a_measured(1,:),'r--',data.time,a_measured(2,:),'b--',data.time,a_measured(3,:),'g--','lineWidth',2.0);
axis tight;
grid on; 

figure();
subplot(211);
plot(data.time,f_fp(:,4),'r',data.time,f_fp(:,5),'b',data.time,f_fp(:,6),'g', 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('ForcePlate Force (N)','FontSize',18);
title('Predicted vs Actual Wrench in forceplate','FontSize',15);
        legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
hold on; 
plot(data.time,f_measured(4,:)','r--',data.time,f_measured(5,:)','b--',data.time,f_measured(6,:),'g--','lineWidth',2.0);
axis tight;
grid on;

subplot(212);
plot(data.time,f_fp(:,1),'r',data.time,f_fp(:,2),'b',data.time,f_fp(:,3),'g', 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('ForcePlate Momment (Nm)','FontSize',18);
legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
hold on;
plot(data.time,f_measured(1,:)','r--',data.time,f_measured(2,:)','b--',data.time,f_measured(3,:),'g--','lineWidth',2.0);
axis tight;
grid on;

 %% Plotting results
save('./experiments/humanFixedBase/resultsFromCheckRNEA.mat', 'tau');