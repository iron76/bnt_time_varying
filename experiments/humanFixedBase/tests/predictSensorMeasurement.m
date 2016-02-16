% predictSensorMeasurement
% script makes a prediction of the sensor measurements of the Human
% Dynamics Estimation Experiment. The sensors are the VICON markers, IMU
% placed on the chest and force place on the bottom of the foot. The method
% loads the sensor measurements (q, dq) and the link sensor transforms 
% from the intermediateDataFiles folder. 
%
% It can be useful exploiting this analysis as a benchmark.
%
% Toolbox requirements:
% - iDynTree - mex
% - Featherstone toolbox (v2)
% 
% Prerequisite : 
% Run synchroniseCaptureData, computeLinkSensorFrames and
% organiseSensorData if running this function for the first time
%
% Author: Naveen Kuppuswamy (naveen.kuppuswamy@iit.it)
% iCub Facility, Istituto Italiano di Tecnologia, 21 January 2016

clear;clc;close all;
%colors = colormap (parula(128));

%% testOptions
% q_i, dq_i and tau_i
plotJointQuantities = true;

% Quatities in link frame associated to sensor measurements (i.e. a2_linear,v2_angular, and f1)
plotLinkQuantities = true; % Link Frame prediction test below will plot same quantities again..

% To analyse XStar_fp_1(t) accuracy (instead of assuming its a constant since joint 1 rotation changes the rotation matrix slightly)
plotTimeVaryingQuantities = false; 

% To run checkRNEA against iDynTree
analyseRNEA = false; % Not yet tested with new refactor (pls retain to false)

% Project sensor quantities in link frame to compare with link quantities
plotLinkFramePrediction = true; 

%% Load Drake model and data
load('./experiments/humanFixedBase/data/humanThreeLinkModelFromURDF_subject1.mat');
load('./experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat','processedSensorData');

t = processedSensorData(1,1).t;
f_capture = processedSensorData(1,1).f_fp_fp;

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

f_measured = processedSensorData(1,1).f_fp_fp(:,tminIndex:tmaxIndex);
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
f_1_1 = zeros(size(q,1),6); % force transmitted from link0 to link1 expressed in link0 frame
a_2_2 = zeros(size(q,1),6); % spatial acceleration link2
v_2_2 = zeros(size(q,1),6); % spatial velocity link2

fprintf('Iterating through time, computing dynamics using RNEA\n');
%humanThreeLink_dmodel.jtype = {'Ry','Ry'};
for i = 1:size(q)
      [tau_i, a_i, v_i, fB_i, f_i] = IDv( humanThreeLink_dmodel, q(i,:), dq(i,:), ddq(i,:));
       tau(i,:) = tau_i';
       a_2_2(i,:) = a_i{2}';
       v_2_2(i,:) = v_i{2}';
       f_1_1(i,:) = f_i{1}';
   
end

if(plotJointQuantities)
    %% Plot raw ID components
    
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    
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

%%  ========================= Sensor prediction ===========================
% Measurements equations are:
% y_2(acc)  = S_lin*(imu_X_2 * a_2_2) + S_ang*(imu_X_2 * v_2_2) x S_lin*(imu_X_2 * v_2_2)
% y_2(gyro) = S_ang*(imu_X_2 * v_2_2)
% y_1(fp)   = fp_XStar_0 * (0_XStar_1* f_1_1 - I_0 * ag)
%
% To predict  measurements we use quantities coming from ID Featherstone


load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');
X_imu_2 = sensorLinkTransforms.X_imu_2;
XStar_fp_0 = sensorLinkTransforms.XStar_fp_0;
XStar_0_1 = sensorLinkTransforms.XStar_0_1;
nT = size(a_2_2,1);

S_lin = [zeros(3) eye(3)];
S_ang = [eye(3) zeros(3)];

 

% IMU PREDICTION
a_imu_imu = S_lin*(X_imu_2*a_2_2') + skew(S_ang*X_imu_2*v_2_2')*(S_lin*X_imu_2*v_2_2'); %skew or cross product are equivalent
omega_imu_imu = S_ang*(X_imu_2*v_2_2');

% FORCE PLATE PREDICTION
a_grav = [0;0;0;0;0;-9.8100]; %Featherstone-like notation

I_c = [0.003 0 0; 0 0.009 0; 0 0 0.012]; %values from URDF file for subject_1
I_0 = createSpatialInertia(I_c,2.057,[0;0;0.026]);

f_0_1 = (XStar_0_1 * f_1_1')';
fg0_0 = repmat( (I_0 * a_grav)',length(f_0_1),1);
f_fp_fp = (XStar_fp_0 * (f_0_1-fg0_0)')'; % -I_0 * g ;

%% Plot link quantities from ID Featherstone
if(plotLinkQuantities)

    
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    
    subplot(211);
    plot1 = plot(data.time,S_lin*a_2_2','lineWidth',2.0); hold on; 
    leg = legend('$a_x$','$a_y$','$a_z$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',15);
    ylabel('Lin Acceleration (m/sec^2)','FontSize',15);
    title('Lin Acceleration and Ang Velocity in frame associated to link2','FontSize',14);
    axis tight; 
    grid on;
    
    subplot(212);
    plot2 = plot(data.time,S_ang*v_2_2','lineWidth',2.0); hold on;
    leg = legend('$w_x$','$w_y$','$w_z$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',15);
    ylabel('Ang Velocity (rad/sec)','FontSize',15);
    axis tight; 
    grid on;
    
    
    figure()
    subplot(211)
    plot3 = plot(data.time,f_1_1(:,4:6),'lineWidth',2.0); hold on;
    leg = legend('$F_x$','$F_y$','$F_z$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',15);
    ylabel('Force (N)','FontSize',15);
    title('Wrench in frame associated to link1 ','FontSize',14);
    axis tight; 
    grid on;
    
    subplot(212)
    plot4 = plot(data.time,f_1_1(:,1:3),'lineWidth',2.0); hold on;
    leg = legend('$\mu_x$','$\mu_y$','$\mu_z$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',15);
    ylabel('Moment (Nm)','FontSize',15);
    axis tight; 
    grid on;    
end

fprintf('Predicting sensor measurement using current sensor-link transforms\n');
    
if(analyseRNEA == true)
    checkRNEA_iDynTree;
end


%% Accelerometer prediction comparison

fig = figure();
% fig = figure('name','xxx');
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');
grid on;


subplot(221)
plot(data.time,a_imu_imu(1,:),data.time,a_imu_imu(2,:),data.time,a_imu_imu(3,:), 'lineWidth',2.0);
leg = legend('$a_x$','$a_y$','$a_z$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',15);
ylabel('Prediction ','FontSize',15);
title('Linear acceleration (m/sec^2)','FontSize',14);
axis tight;
grid on; 

subplot(223)
plot(data.time,a_measured(1,:),data.time,a_measured(2,:),data.time,a_measured(3,:),'lineWidth',2.0);
leg = legend('$a_x$','$a_y$','$a_z$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',15);
ylabel('Actual','FontSize',15);
axis tight;
grid on; 


%% Gyroscope prediction comparison

subplot(222)
plot(data.time,omega_imu_imu(1,:),data.time,omega_imu_imu(2,:),data.time,omega_imu_imu(3,:), 'lineWidth',2.0);
leg = legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',15);
ylabel('Prediction','FontSize',15);
title('Angular velocity (rad/sec)','FontSize',14);
axis tight;
grid on; 

subplot(224)
plot(data.time,omega_measured(1,:),data.time,omega_measured(2,:),data.time,omega_measured(3,:),'lineWidth',2.0);
leg = legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',15);
ylabel('Actual','FontSize',15);
axis tight;
grid on; 

axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
text(0.5, 0.99,'\bf Prediction comparison (IMU frame)','HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',14);

%% Force prediction comparison 

figure();
subplot(221);
plot(data.time,f_fp_fp(:,4),data.time,f_fp_fp(:,5),data.time,f_fp_fp(:,6), 'lineWidth',2.0);
leg = legend('$F_x$','$F_y$','$F_z$','Location','southeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',15);
ylabel('Prediction','FontSize',15);
title('Force (N)','FontSize',15);
hold on; 
axis tight;
grid on;

subplot(223)
plot(data.time,f_measured(4,:)',data.time,f_measured(5,:)',data.time,f_measured(6,:),'lineWidth',2.0);
leg = legend('$F_x$','$F_y$','$F_z$','Location','southeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',15);
ylabel('Actual','FontSize',15);
axis tight;
grid on;

%% Moment prediction comparison 

subplot(222);
plot(data.time,f_fp_fp(:,1),data.time,f_fp_fp(:,2),data.time,f_fp_fp(:,3), 'lineWidth',2.0);
leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',15);
ylabel('Prediction','FontSize',15);
title('Moment (Nm)','FontSize',15);
hold on;
axis tight;
grid on;

subplot(224)
plot(data.time,f_measured(1,:)',data.time,f_measured(2,:)',data.time,f_measured(3,:),'lineWidth',2.0);
leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',15);
ylabel('Actual','FontSize',15);
axis tight;
grid on;


axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
text(0.5, 0.99,'\bf Prediction comparison (Fp frame)','HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',14);


