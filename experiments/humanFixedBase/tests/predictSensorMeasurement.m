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

%% testOptions
% q_i, qDot_i and tau_i
plotJointQuantities = true;

% Quatities in link frame associated to sensor measurements (i.e. a2_linear,v2_angular, and f1)
plotLinkQuantities = true; 

% To analyse X_fp_1(t) accuracy (instead of assuming its a constant since joint 1 rotation changes the rotation matrix slightly)
plotTimeVaryingQuantities = false; 

% To run checkRNEA against iDynTree
analyseRNEA = false; 

% Project sensor quantities in link frame to compare with link quantities
plotLinkFramePrediction = true; 

%% Load Drake model
load('humanThreeLinkModelFromURDF_subject1.mat');

%% Load acquired data
load('./experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat','processedSensorData');

t = processedSensorData(1,1).t;
f_capture = processedSensorData(1,1).f_fp;

% tminIndex = 1;
% tmaxIndex = length(t);
% data.min_time = t(1);
% data.max_time = t(end);
% data.nsamples = length(t);
% data.time = t;

[~,chosenF_ID] = max(f_capture(:,1));
[~,tminIndex] = max(f_capture(chosenF_ID,1:round(end/2)));
[~,tmaxIndex] = max(f_capture(chosenF_ID,round(end/2):end));
totPointsInConsideration = tmaxIndex - tminIndex;
tminIndex = tminIndex + round(0.1*totPointsInConsideration);
tmaxIndex = tmaxIndex + round(length(t)/2) - round(0.1*totPointsInConsideration);

data.min_time = t(tminIndex);
data.max_time = t(tmaxIndex);
data.nsamples = tmaxIndex - tminIndex;
data.time = t(tminIndex:tmaxIndex);

f_measured = processedSensorData(1,1).f_fp(:,tminIndex:tmaxIndex);
a_measured = processedSensorData(1,1).a_imu(:,tminIndex:tmaxIndex);
omega_measured = processedSensorData(1,1).omega_imu(:,tminIndex:tmaxIndex);

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

%GET JOINT ACCELERATIONS 
ddq1 = processedSensorData.ddq1;
ddq1 = ddq1(tminIndex:tmaxIndex,:); %window filter
ddq2 = processedSensorData.ddq2;
ddq2 = ddq2(tminIndex:tmaxIndex,:); %window filter

q  = [q1,q2];
dq  = [dq1,dq2];
ddq  = [ddq1,ddq2];

%% Computing tau using Newton-Euler with Featherstone toolbox

tau = zeros(size (q));
f_1 = zeros(size(q,1),6);
a_2 = zeros(size(q,1),6);
v_2 = zeros(size(q,1),6);

fprintf('Iterating through time, computing dynamics using RNEA\n');
%humanThreeLink_dmodel.jtype = {'Ry','Ry'};
for i = 1:size(q)
      [tau_i, a_i, v_i, fB_i, f_i] = ID( humanThreeLink_dmodel, q(i,:), dq(i,:), ddq(i,:));
      tau(i,:) = tau_i';
       a_2(i,:) = a_i{2}';
       v_2(i,:) = v_i{2}';
       f_1(i,:) = f_i{1}';
      % f(i,:) = f_i';   
end

if(plotJointQuantities)
    %% raw ID components
    figure();
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
end

%% Computing the sensor measurement in its frame

load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');
X_imu_2 = sensorLinkTransforms.X_imu_2;
XStar_fp_1ini = sensorLinkTransforms.XStar_fp_1ini;
XStar_fp_1t = sensorLinkTransforms.XStar_fp_1t;
nT = size(a_2,1);

X_2_imu = sensorLinkTransforms.X_2_imu;
XStar_1_fpini = sensorLinkTransforms.XStar_1_fpini;
XStar_1_fpt = sensorLinkTransforms.XStar_1_fpt;

S_lin = [zeros(3) eye(3)];
S_ang = [eye(3) zeros(3)];

%% Computed IMU Prediction
a_imu = S_lin*(X_imu_2*a_2') + skew(S_ang*X_imu_2*v_2')*(S_lin*X_imu_2*v_2');
omega_imu = S_ang*(X_imu_2*v_2');

%% Computed Forceplate Prediction
%f_fp = (XStar_0_fp *fB_1')';
f_fp = (XStar_fp_1ini *f_1')';

%% Computing inverse predictions (in link frames)
a2_pred = X_2_imu * (pinv(S_lin)*a_measured);
v2_pred = X_2_imu * (pinv(S_ang)*omega_measured);
f1_pred = (XStar_1_fpini * f_measured)';


if(plotLinkQuantities)
    figure();
    plot1 = plot(data.time,S_lin*a_2','lineWidth',2.0); hold on;
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Joint 2 acceleration (m/sec^2)','FontSize',18);
    title('Joint 2 acceleration');
    legend('$a_x$','$a_y$','$az$','Location','northeast');
    axis tight; grid on;
    
    figure();
    plot1 = plot(data.time,S_ang*v_2','lineWidth',2.0); hold on;
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Joint 2 velocity (rad/sec)','FontSize',18);
    title('Joint 2 velocity');
    legend('$w_x$','$w_y$','$wz$','Location','northeast');
    axis tight; grid on;
    
    figure()
    subplot(211)
    plot1 = plot(data.time,f_1(:,4:6),'lineWidth',2.0); hold on;
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Wrench Force(N)','FontSize',18);
    title('Link1 transmitted wrench to base');
    legend('$F_x$','$F_y$','$Fz$','Location','northeast');
    axis tight; grid on;
    subplot(212)
    plot1 = plot(data.time,f_1(:,1:3),'lineWidth',2.0); hold on;
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Wrench Moment(Nm)','FontSize',18);
    legend('$\mu_x$','$\mu_y$','$\mu_z$','Location','northeast');
    axis tight; grid on;    
end

fprintf('Predicting sensor measurement using current sensor-link transforms\n');
    
if(analyseRNEA == true)
    checkRNEA_iDynTree;
end


%% Accelerometer prediction comparison
figure();
subplot(211)
plot(data.time,a_imu(1,:),'r',data.time,a_imu(2,:),'b',data.time,a_imu(3,:),'g', 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Accelerometer prediction (m/sec^2)','FontSize',18);
title('Predicted vs Actual Accelerometer measurements');
legend('$a_x$','$a_y$','$a_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on; 
%hold on;
subplot(212)
plot(data.time,a_measured(1,:),'r',data.time,a_measured(2,:),'b',data.time,a_measured(3,:),'g','lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Accelerometer measurement (m/sec^2)','FontSize',18);
legend('$a_x$','$a_y$','$a_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on; 

%% Gyroscope prediction comparison
%omega
figure();
subplot(211)
plot(data.time,omega_imu(1,:),'r',data.time,omega_imu(2,:),'b',data.time,omega_imu(3,:),'g', 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Gyroscope prediction (rad/sec)','FontSize',18);
title('Predicted vs Actual Gyroscope measurements');
legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on; 
%hold on;
subplot(212)
plot(data.time,omega_measured(1,:),'r',data.time,omega_measured(2,:),'b',data.time,omega_measured(3,:),'g','lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Gyroscope measurement (rad/sec)','FontSize',18);
title('Prediction vs Actual Gyroscope measurements');
legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on; 

%% Forceplate prediction comparison (assuming constant XStar_1_fp)
figure();
subplot(211);
plot(data.time,f_fp(:,4),'r',data.time,f_fp(:,5),'b',data.time,f_fp(:,6),'g', 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('ForcePlate Force (N)','FontSize',18);
title('Predicted vs Actual Force in forceplate','FontSize',15);
        legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
hold on; 
axis tight;
grid on;

subplot(212)
plot(data.time,f_measured(4,:)','r',data.time,f_measured(5,:)','b',data.time,f_measured(6,:),'g','lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('ForcePlate Force (N)','FontSize',18);
        legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on;


figure()
subplot(211);
plot(data.time,f_fp(:,1),'r',data.time,f_fp(:,2),'b',data.time,f_fp(:,3),'g', 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('ForcePlate Momment (Nm)','FontSize',18);
legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
title('Predicted vs Actual Moment in forceplate','FontSize',15);
hold on;
axis tight;
grid on;

subplot(212)
plot(data.time,f_measured(1,:)','r',data.time,f_measured(2,:)','b',data.time,f_measured(3,:),'g','lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('ForcePlate Momment (Nm)','FontSize',18);
legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on;


%% Plotting the predicted quantities in a time varying case (Because XStar_1_fp is not exactly a constant)
if(plotTimeVaryingQuantities)
    
    for i = 1:length(data.time)
        f_fpt(i,:) = (XStar_fp_1t{i} *f_1(i,:)')';
    end
    figure();
    subplot(211);
    plot(data.time,f_fpt(:,4),'r',data.time,f_fpt(:,5),'b',data.time,f_fpt(:,6),'g', 'lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('ForcePlate Force (N)','FontSize',18);
    title('Predicted vs Actual Force in forceplate (timeVarying)','FontSize',15);
            legend('$F_x$','$F_y$','$F_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    hold on; 
    axis tight;
    grid on;

    subplot(212)
    plot(data.time,f_measured(4,:)','r',data.time,f_measured(5,:)','b',data.time,f_measured(6,:),'g','lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('ForcePlate Force (N)','FontSize',18);
            legend('$F_x$','$F_y$','$F_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    axis tight;
    grid on;


    figure()
    subplot(211);
    plot(data.time,f_fpt(:,1),'r',data.time,f_fpt(:,2),'b',data.time,f_fpt(:,3),'g', 'lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('ForcePlate Momment (Nm)','FontSize',18);
    legend('$M_x$','$M_y$','$M_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    title('Predicted vs Actual Moment in forceplate (timeVarying)','FontSize',15);
    hold on;
    axis tight;
    grid on;

    subplot(212)
    plot(data.time,f_measured(1,:)','r',data.time,f_measured(2,:)','b',data.time,f_measured(3,:),'g','lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('ForcePlate Momment (Nm)','FontSize',18);
    legend('$M_x$','$M_y$','$M_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    axis tight;
    grid on;
    
end

%% Plotting the comparisons between link frame predictions and link frame quantities
if(plotLinkFramePrediction)
    figure();
    subplot(211)
    plot(data.time,a2_pred(4,:),'r',data.time,a2_pred(5,:),'b',data.time,a2_pred(6,:),'g', 'lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Link LinAcceleration prediction (m/sec^2)','FontSize',18);
    title('Predicted vs Actual Link Frame LinAcceleration');
    legend('$a_x$','$a_y$','$a_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    axis tight;
    grid on; 
    %hold on;
    subplot(212)
    plot(data.time,a_2(:,4),'r',data.time,a_2(:,5),'b',data.time,a_2(:,6),'g','lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Link LinAcceleration actual (m/sec^2)','FontSize',18);
    legend('$a_x$','$a_y$','$a_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    axis tight;
    grid on; 

    %omega
    figure();
    subplot(211)
    plot(data.time,v2_pred(1,:),'r',data.time,v2_pred(2,:),'b',data.time,v2_pred(3,:),'g', 'lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Link AngVelocity prediction (rad/sec)','FontSize',18);
    title('Predicted vs Actual Link Frame AngVelocity');
    legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    axis tight;
    grid on; 
    %hold on;
    subplot(212)
    plot(data.time,v_2(:,1),'r',data.time,v_2(:,2),'b',data.time,v_2(:,3),'g','lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Link AngVelocity Actual (rad/sec)','FontSize',18);
    legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    axis tight;
    grid on; 

    figure();
    subplot(211);
    plot(data.time,f1_pred(:,4),'r',data.time,f1_pred(:,5),'b',data.time,f1_pred(:,6),'g', 'lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Link Transmitted Force (N)','FontSize',18);
    title('Predicted vs Actual Force transmitted to base by link1','FontSize',15);
            legend('$F_x$','$F_y$','$F_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    hold on; 
    axis tight;
    grid on;

    subplot(212)
    plot(data.time,f_1(:,4),'r',data.time,f_1(:,5),'b',data.time,f_1(:,6),'g','lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('ForcePlate Force (N)','FontSize',18);
            legend('$F_x$','$F_y$','$F_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    axis tight;
    grid on;


    figure()
    subplot(211);
    plot(data.time,f1_pred(:,1),'r',data.time,f1_pred(:,2),'b',data.time,f1_pred(:,3),'g', 'lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('ForcePlate Momment (Nm)','FontSize',18);
    legend('$M_x$','$M_y$','$M_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    title('Predicted vs Actual Moment  transmitted to base by link1','FontSize',15);
    hold on;
    axis tight;
    grid on;

    subplot(212)
    plot(data.time,f_1(:,1)','r',data.time,f_1(:,2)','b',data.time,f_1(:,3),'g','lineWidth',2.0);
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('ForcePlate Momment (Nm)','FontSize',18);
    legend('$M_x$','$M_y$','$M_z$','Location','northeast');
            set(legend,'Interpreter','latex');
            set(legend,'FontSize',20);
    axis tight;
    grid on;
end
