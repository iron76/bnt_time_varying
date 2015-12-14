

% Check with iDynTree lib

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

leg_pos_foot = leg_RelTransf_foot.getPosition().toMatlab(); %matlab object
leg_R_foot = leg_RelTransf_foot.getRotation().toMatlab(); %matlab object
leg_X_foot = leg_RelTransf_foot.asAdjointTransform().toMatlab(); %matlab object


% transformation from leg to torso
torso_RelTransf_leg = model.getRelativeTransform('torso','leg'); %c++ object

torso_pos_leg  = torso_RelTransf_leg.getPosition().toMatlab(); %matlab object
torso_R_leg = torso_RelTransf_leg.getRotation().toMatlab(); %matlab object
torso_X_leg = torso_RelTransf_leg.asAdjointTransform().toMatlab(); %matlab object


%% Plot tau using Newton-Euler with iDynTree lib

tauIdyntree = iDynTreeID(model,q,dq,ddq);

% fig = figure();
% axes1 = axes('Parent',fig,'FontSize',16);
% box(axes1,'on');
% hold(axes1,'on');
% grid on;
% 
% plot1 = plot(data.time,tauIdyntree(:,1),'lineWidth',3.0, 'Parent',axes1);hold on;
% set(plot1,'color',[1 0 0]);
% plot2 = plot(data.time,tauIdyntree(:,2),'lineWidth',3.0, 'Parent',axes1);hold on;
% set(plot2,'color',[0 0.498039215803146 0]);
% leg = legend('$\tau_1$','$\tau_2$');
% set(leg,'Interpreter','latex');
% set(leg,'FontSize',18);
% xlabel('Time [s]','FontSize',20);
% ylabel('\tau_{IdynTree}[Nm]','FontSize',20);
% axis tight; hold on;   

%%  Comparing error

error = norm(tau-tauIdyntree);
fprintf('Inverse dynamics compared with an error of %d', error);

 % Comparing plots Featherstone/IdynTree

figure;
subplot(511);
plot1 = plot(data.time,(180/pi)*q1,'lineWidth',2.0); hold on;
set(plot1,'color',[1 0 0]);
plot2= plot(data.time,(180/pi)*q2,'lineWidth',2.0); hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$q_1$','$q_2$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',12);
ylabel('Joint Ang [deg]','FontSize',12);
axis tight;
grid on;
 
subplot(512);
plot1 = plot(data.time,(180/pi)*dq1,'lineWidth',2.0); hold on;
set(plot1,'color',[1 0 0]);
plot2= plot(data.time,(180/pi)*dq2,'lineWidth',2.0); hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\dot q_1$','$\dot q_2$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',12);
ylabel('Joint Vel [deg/s]','FontSize',12);
axis tight;
grid on;

subplot(513);
plot1 = plot(data.time,(180/pi)*ddq1,'lineWidth',2.0);hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,(180/pi)*ddq2,'lineWidth',2.0);hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\ddot q_1$','$\ddot q_2$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',12);
ylabel('Joint Acc [deg/s^2]','FontSize',12);
axis tight; 
grid on;

subplot(514);
plot1 = plot(data.time,tau(:,1), 'lineWidth',2.0); hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,tau(:,2), 'lineWidth',2.0); hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\tau_1$','$\tau_2$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',12);
ylabel('\tau_{RNEA,Feath.} [Nm]','FontSize',12);
axis tight;
grid on;

subplot(515);
plot1 = plot(data.time,tauIdyntree(:,1), 'lineWidth',2.0); hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,tauIdyntree(:,2), 'lineWidth',2.0); hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\tau_1$','$\tau_2$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',12);
ylabel('\tau_{IdynTree}[Nm]','FontSize',12);
axis tight;
grid on;
