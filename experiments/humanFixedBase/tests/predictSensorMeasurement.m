% predictSensorMeasurement
% function makes a prediction of the sensor measurements of the Human
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
plotLinkQuantities = false; % Link Frame prediction test below will plot same quantities again..

% To analyse X_fp_1(t) accuracy (instead of assuming its a constant since joint 1 rotation changes the rotation matrix slightly)
plotTimeVaryingQuantities = false; 

% To run checkRNEA against iDynTree
analyseRNEA = false; % Not yet tested with new refactor (pls retain to false)

% Project sensor quantities in link frame to compare with link quantities
plotLinkFramePrediction = true; 

%% Load Drake model
load('./experiments/humanFixedBase/data/humanThreeLinkModelFromURDF_subject1.mat');

%% Load acquired data
load('./experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat','processedSensorData');

%% Load Link Sensor Transforms
load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');

X_imu_2 = sensorLinkTransforms.X_imu_2;
XStar_fp_0 = sensorLinkTransforms.XStar_fp_0;
XStar_0_1 = sensorLinkTransforms.XStar_0_1;

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

tau = zeros(size (q)); % joint torques
f_1 = zeros(size(q,1),6); % force transmitted to link 0 expressed in link0 frame
a_2 = zeros(size(q,1),6); % spatial acceleration link 2
v_2 = zeros(size(q,1),6); % spatial velocity link 2

f_0_1 = zeros(size(q,1),6);

fprintf('Iterating through time, computing dynamics using RNEA\n');
%humanThreeLink_dmodel.jtype = {'Ry','Ry'};
for i = 1:size(q)
      [tau_i, a_i, v_i, fB_i, f_i] = IDv( humanThreeLink_dmodel, q(i,:), dq(i,:), ddq(i,:));
      tau(i,:) = tau_i';
       a_2(i,:) = a_i{2}';
       v_2(i,:) = v_i{2}';
       f_1(i,:) = f_i{1}';
       
       if(iscell(XStar_0_1)==1)
           f_0_1(i,:) = ( XStar_0_1{i}* f_1(i,:)')';
       end
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
    title('RNEA inputs : q, \dot{q}, \ddot{q} (Actual)');

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

nT = size(a_2,1);

S_lin = [zeros(3) eye(3)];
S_ang = [eye(3) zeros(3)];

%% Computed IMU Prediction
a_imu_s = S_lin*(X_imu_2*a_2') + skew(S_ang*X_imu_2*v_2')*(S_lin*X_imu_2*v_2');
a_imu = S_lin*(X_imu_2*a_2') + cross((S_ang*X_imu_2*v_2'),(S_lin*X_imu_2*v_2'));

omega_imu = S_ang*(X_imu_2*v_2');

%% Computed Forceplate Prediction


a_grav = [0;0;0;0;0; -9.8100]; %Featherstone-like notation

I_c = [0.003 0 0; 0 0.009 0; 0 0 0.012]; %values from URDF file for subject_1
I_0 = createSpatialInertia(I_c,2.057,[0;0;0.026]);


fg0_0 = repmat( (I_0*a_grav)',length(f_0_1),1);
f_fp = (XStar_fp_0 * (f_0_1-fg0_0)')'; %-I0 g ;

if(plotLinkQuantities)
    figure();
    subplot(211);
    plot1 = plot(data.time,S_lin*a_2','lineWidth',2.0); hold on;
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Linacceleration (m/sec^2)','FontSize',18);
    title('Link 2 LinAcceleration and AngVelocity');
    legend('$a_x$','$a_y$','$az$','Location','northeast');
    axis tight; grid on;
    subplot(212);
    plot1 = plot(data.time,S_ang*v_2','lineWidth',2.0); hold on;
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Angvelocity (rad/sec)','FontSize',18);
    legend('$w_x$','$w_y$','$wz$','Location','northeast');
    axis tight; grid on;
    
    figure()
    subplot(211)
    plot1 = plot(data.time,f_1(:,4:6),'lineWidth',2.0); hold on;
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Wrench Force(N)','FontSize',18);
    title('Link1 transmitted wrench');
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
subplot(221)
plot(data.time,a_imu(1,:),'r',data.time,a_imu(2,:),'b',data.time,a_imu(3,:),'g', 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('LinAcceleration prediction(m/sec^2)','FontSize',18);
title('Accelerometer');
legend('$a_x$','$a_y$','$a_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on; 
%hold on;
subplot(223)
plot(data.time,a_measured(1,:),'r',data.time,a_measured(2,:),'b',data.time,a_measured(3,:),'g','lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('LinAccleration Actual(m/sec^2)','FontSize',18);
legend('$a_x$','$a_y$','$a_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on; 

%% Gyroscope prediction comparison
%omega
subplot(222)
plot(data.time,omega_imu(1,:),'r',data.time,omega_imu(2,:),'b',data.time,omega_imu(3,:),'g', 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('AngVelocity Prediction(rad/sec)','FontSize',18);
title('Gyroscope');
legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on; 
%hold on;
subplot(224)
plot(data.time,omega_measured(1,:),'r',data.time,omega_measured(2,:),'b',data.time,omega_measured(3,:),'g','lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('AngVelocity Actual(rad/sec)','FontSize',18);
legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on; 

%% Forceplate prediction comparison (assuming constant XStar_1_fp)
figure();
subplot(221);
plot(data.time,f_fp(:,4),'r',data.time,f_fp(:,5),'b',data.time,f_fp(:,6),'g', 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Force Prediction(N)','FontSize',18);
title('Force in forceplate','FontSize',15);
        legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
hold on; 
axis tight;
grid on;

subplot(223)
plot(data.time,f_measured(4,:)','r',data.time,f_measured(5,:)','b',data.time,f_measured(6,:),'g','lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Force Actual(N)','FontSize',18);
        legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on;


%figure()
subplot(222);
plot(data.time,f_fp(:,1),'r',data.time,f_fp(:,2),'b',data.time,f_fp(:,3),'g', 'lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Momment Prediction(Nm)','FontSize',18);
legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
title('Moment in forceplate','FontSize',15);
hold on;
axis tight;
grid on;

subplot(224)
plot(data.time,f_measured(1,:)','r',data.time,f_measured(2,:)','b',data.time,f_measured(3,:),'g','lineWidth',2.0);
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',18);
ylabel('Momment Actual(Nm)','FontSize',18);
legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
axis tight;
grid on;




