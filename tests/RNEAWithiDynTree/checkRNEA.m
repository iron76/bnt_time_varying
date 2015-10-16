
%% CHECK RNEA 
% provides the RNEA analysis assuming no external forces.
%It can be useful exploiting this analysis as a benchmark.
% Reminder: if Drake path (MATLAB/Drake-Distro/drake/ThirdParty/Spatia)
% is not included in Matlab path, the recalled ID function  is that one of
% Featherstone, In order to avoid issue, be sure using Featherstone ID 
% function!
% 
%%Toolbox requirements:
%% - iDynTree - mex
%% - Featherstone toolbox (v2) with ID corrected as in bnt_time_varying repository
%%
clear;clc;


%% Load Drake model
load('humanThreeLinkModelFromURDF_subject1.mat');
% humanThreeLink_dmodel.gravity = [0;0;-9.81]; 
% humanThreeLink_dmodel.f_ext_map_from = [];
% humanThreeLink_dmodel.f_ext_map_to = [];

%% Check model imported from Drake

leg_R_foot = humanThreeLink_dmodel.Xtree{1}(1:3,1:3)
torso_R_leg = humanThreeLink_dmodel.Xtree{2}(1:3,1:3);

leg_r_foot = skew(-leg_R_foot' * humanThreeLink_dmodel.Xtree{1}(4:6,1:3))


%% Load acquired data
load('preProcessedSensorData.mat','processedSensorData');

t = processedSensorData(1,1).t;
f_temp = processedSensorData(1,1).ftx;
[val,chosenF_ID] = max(f_temp(:,1));
[val,tminIndex] = max(f_temp(chosenF_ID,1:round(end/2)));
[val,tmaxIndex] = max(f_temp(chosenF_ID,round(end/2):end));
totPointsInConsideration = tmaxIndex - tminIndex;
tminIndex = tminIndex + round(0.1*totPointsInConsideration);
tmaxIndex = tmaxIndex + round(length(t)/2) - round(0.1*totPointsInConsideration);
data.min_time = t(tminIndex);
data.max_time = t(tmaxIndex);
data.nsamples = tmaxIndex - tminIndex;
data.time = t(tminIndex:tmaxIndex);

%% Getting variables

%GET JOINT ANGLES
q1 = processedSensorData.q1;
q1 = q1(2000:18606,:); %window filter
q2 = processedSensorData.q2;
q2 = q2(2000:18606,:); %window filter


%GET JOINT VELOCITIES 
% dq1 = processedSensorData.dq1;
% dq1 = dq1(2000:18606,:); %window filter
% dq2 = processedSensorData.dq2;
% dq2 = dq2(2000:18606,:); %window filter
dq1 = zeros (length(q1),1);
dq2 = zeros (length(q2),1);


%GET JOINT ACCELERATIONS
% % ddq1 = processedSensorData.ddq1;
% % ddq1 = ddq1(2000:18606,:); %window filter
% % ddq2 = processedSensorData.ddq2;
% % ddq2 = ddq2(2000:18606,:); %window filter
ddq1 = zeros (length(q1),1);
ddq2 = zeros (length(q2),1);

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

%% Plot tau using Newton-Euler with Featherstone toolbox

q  = [q1,q2];
dq  = [dq1,dq2];
ddq  = [ddq1,ddq2];

tau = zeros(size (q));
% a = cell (size(q));
% fB_i = cell (size(q));
% f_i = cell (size(q));

[tau_i, a_i, fB_i, f_i] = ID( humanThreeLink_dmodel, [45/180*pi;0], [0;0], [0;0]);

for i = 1:size(q)
      [tau_i, a_i, fB_i, f_i] = ID( humanThreeLink_dmodel, q(i,:), dq(i,:), ddq(i,:));
      tau(i,:) = tau_i';
%       a(i,:) = a_i';
%       fB(i,:) = fB_i';
%       f(i,:) = f_i';   
end

fig = figure();
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');
grid on;

plot1 = plot(data.time,tau(:,1),'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,tau(:,2),'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\tau_1$','$\tau_2$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Torque_{RNEA-TDA} [Nm]','FontSize',20);
axis tight; hold on;   


%% Check with iDynTree

model = iDynTree.DynamicsComputations();

if ~model.loadRobotModelFromFile('threeLinkHuman_subject1.urdf')
end

dofName = model.getDescriptionOfDegreesOfFreedom(); %c++ object
NumLink =model.getNrOfLinks(); %c++ object
NumFrame = model.getNrOfFrames(); %c++ object

for i = 0:(NumFrame-1)
    model.getFrameName(i);
end


for i = 0:(NumLink-1)
    inertia =  model.getLinkInertia(i); %c++ object
    linkMass = inertia.getMass(); %c++ object
    linkCOM = inertia.getCenterOfMass().toMatlab();
end


% transformation from foot to leg
leg_RelTransf_foot = model.getRelativeTransform('leg','foot'); %c++ object

leg_pos_foot = leg_RelTransf_foot.getPosition().toMatlab() %matlab object
leg_R_foot = leg_RelTransf_foot.getRotation().toMatlab() %matlab object
leg_X_foot = leg_RelTransf_foot.asAdjointTransform().toMatlab(); %matlab object


% transformation from leg to torso
torso_RelTransf_leg = model.getRelativeTransform('torso','leg'); %c++ object

torso_pos_leg  = torso_RelTransf_leg.getPosition().toMatlab(); %matlab object
torso_R_leg = torso_RelTransf_leg.getRotation().toMatlab(); %matlab object
torso_X_leg = torso_RelTransf_leg.asAdjointTransform().toMatlab(); %matlab object

%% Plot tau using Newton-Euler with iDynTree lib

tauIdyntree = iDynTreeID(model,q,dq,ddq);

fig = figure();
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');
grid on;

plot1 = plot(data.time,tauIdyntree(:,1),'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,tauIdyntree(:,2),'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\tau_1$','$\tau_2$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Torque_{IdynTree-TDA}[Nm]','FontSize',20);
axis tight; hold on;   


%%  COMPARING ERROR

error = tau + tauIdyntree; 
% sign to be fixed!

fig = figure();
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');
grid on;

plot1 = plot(data.time,error(:,1),'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,error(:,2),'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\tau_1$','$\tau_2$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Torque_{IdynTree-TDA}[Nm]','FontSize',20);
axis tight; hold on;   


NormError = norm(tau+tauIdyntree);


 %% COMPARING ANGLES VELOCITY AND TORQUES_RNEA
% % 
% % figure;
% % subplot(311);
% % plot1 = plot(data.time,(180/pi)*q1,'lineWidth',2.0); hold on;
% % set(plot1,'color',[1 0 0]);
% % plot2= plot(data.time,(180/pi)*q2,'lineWidth',2.0); hold on;
% % set(plot2,'color',[0 0.498039215803146 0]);
% % leg = legend('$q_1$','$q_2$','Location','northeast');
% % set(leg,'Interpreter','latex');
% % set(leg,'FontSize',18);
% % xlabel('Time [s]','FontSize',15);
% % ylabel('Joint Angle [deg]','FontSize',15);
% % axis tight;
% % grid on;
% %  
% % subplot(312);
% % plot1 = plot(data.time,(180/pi)*dq1,'lineWidth',2.0); hold on;
% % set(plot1,'color',[1 0 0]);
% % plot2= plot(data.time,(180/pi)*dq2,'lineWidth',2.0); hold on;
% % set(plot2,'color',[0 0.498039215803146 0]);
% % leg = legend('$\dot q_1$','$\dot q_2$','Location','northeast');
% % set(leg,'Interpreter','latex');
% % set(leg,'FontSize',18);
% % xlabel('Time [s]','FontSize',15);
% % ylabel('Joint Velocity [deg/s]','FontSize',15);
% % axis tight;
% % grid on;
% % 
% % % subplot(413);
% % % plot1 = plot(data.time,(180/pi)*ddq1,'lineWidth',2.0); hold on;
% % % set(plot1,'color',[1 0 0]);
% % % plot2= plot(data.time,(180/pi)*ddq2,'lineWidth',2.0); hold on;
% % % set(plot2,'color',[0 0.498039215803146 0]);
% % % leg = legend('$\ddot q_1$','$\ddot q_2$','Location','northeast');
% % % set(leg,'Interpreter','latex');
% % % set(leg,'FontSize',18);
% % % xlabel('Time [s]','FontSize',15);
% % % ylabel('Joint acc [deg/s^2]','FontSize',15);
% % % axis tight;
% % % grid on;
% %  
% % subplot(313);
% % plot1 = plot(data.time,tau(:,1), 'lineWidth',2.0); hold on;
% % set(plot1,'color',[1 0 0]);
% % plot2 = plot(data.time,tau(:,2), 'lineWidth',2.0); hold on;
% % set(plot2,'color',[0 0.498039215803146 0]);
% % leg = legend('$\tau_1$','$\tau_2$','Location','northeast');
% % set(leg,'Interpreter','latex');
% % set(leg,'FontSize',18);
% % xlabel('Time [s]','FontSize',15);
% % ylabel('Torque_{RNEA-TDA} [Nm]','FontSize',15);
% % axis tight;
% % grid on;
