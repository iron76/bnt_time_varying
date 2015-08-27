close all;

subjectID = 3;
trialID = 1;

%% load dataset
load(sprintf('./experiments/humanFixedBase/savedBERDYresult_subj%d_trial%d.mat',subjectID,trialID),'res','data','myPNEA');

%% plot joint angles
fig = figure(1);
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');

plot1 = plot(data.time,(180/pi)*data.q1-180,'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,(180/pi)*data.q2-180,'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$q_1$','$q_2$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Joint Angle [deg]','FontSize',20);
axis tight; 

%% plot joint velocities
fig = figure(2);
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');

plot1 = plot(data.time,(180/pi)*data.dq1-180,'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,(180/pi)*data.dq2-180,'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\dot q_1$','$\dot q_2$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Joint Velocity [deg/s]','FontSize',20);
axis tight; 

%% plot IMU accelerometer
[RNEA_imu,RNEA_ftx,RNEA_t] = sensorDataConsistencyCheck(subjectID,trialID,'noplots');

fig = figure(3);
%subplot(2,1,1);
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');

plot1 = plot(data.time,data.y_imu(1,:)','lineWidth',3.0, 'Parent',axes1);hold on;
set(plot1,'color',[0 0.447058826684952 0.74117648601532]);
plot2 = plot(data.time,data.y_imu(2,:)','lineWidth',3.0, 'Parent',axes1);hold on;
set(plot2,'color',[0.749019622802734 0 0.749019622802734]);
plot3 = plot(data.time,data.y_imu(3,:)','lineWidth',3.0, 'Parent',axes1);hold on;
set(plot3,'color',[0.929411768913269 0.694117665290833 0.125490203499794]);
h = plot(RNEA_t,RNEA_imu(:,1:3),'--','lineWidth',2.0);
set(h,'color',[0 0 0]);
leg = legend('$a_x$','$a_y$','$a_z$','RNEA');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Acceleration [m/s^2]','FontSize',20);
axis tight; 
%% plot IMU gyro

fig = figure(4);
%subplot(2,1,2);
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');

plot1 = plot(data.time,data.y_imu(4,:)','lineWidth',2.5,'Parent',axes1);hold on;
set(plot1,'color',[0 0.447058826684952 0.74117648601532]);
plot2 = plot(data.time,data.y_imu(5,:)','lineWidth',2.5,'Parent',axes1);hold on;
set(plot2,'color',[0.749019622802734 0 0.749019622802734]);
plot3 = plot(data.time,data.y_imu(6,:)','lineWidth',2.5,'Parent',axes1);hold on;
set(plot3,'color',[0.929411768913269 0.694117665290833 0.125490203499794]);
h = plot(RNEA_t,RNEA_imu(:,4:6),'--','lineWidth',2.0);
set(h,'color',[0 0 0]);
leg = legend('$\omega_x$','$\omega_y$','$\omega_z$','RNEA');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Angular Velocity [rad/s]','FontSize',20);
axis tight; 

%% plot forcePlate force

fig = figure(5);
%subplot(2,1,1);
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');

%plot(data.time,data.y_ftx(1:3,:)','lineWidth',1.0);hold on;
plot1 = plot(data.time,data.y_ftx(1,:)','lineWidth',3.0,'Parent',axes1);hold on;
set(plot1,'color',[0 0.447058826684952 0.74117648601532]);
plot2 = plot(data.time,data.y_ftx(2,:)','lineWidth',3.0,'Parent',axes1);hold on;
set(plot2,'color',[0.749019622802734 0 0.749019622802734]);
plot3 = plot(data.time,data.y_ftx(3,:)','lineWidth',3.0,'Parent',axes1);hold on;
set(plot3,'color',[0.929411768913269 0.694117665290833 0.125490203499794]);

%plot(RNEA_t,RNEA_ftx(:,1:3),'--');
%set(h,'color',[0.5 0.5 0.5]);
leg = legend('$f_x$','$f_y$','$f_z$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Force [N]','FontSize',20);
axis tight; 

%% plot forcePlate torque

fig = figure(6);
%subplot(2,1,2);
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');

%plot(data.time,data.y_ftx(4:6,:)','lineWidth',1.0);hold on;
plot1 = plot(data.time,data.y_ftx(4,:)','lineWidth',3.0,'Parent',axes1);hold on;
set(plot1,'color',[0 0.447058826684952 0.74117648601532]);
plot2 = plot(data.time,data.y_ftx(5,:)','lineWidth',3.0,'Parent',axes1);hold on;
set(plot2,'color',[0.749019622802734 0 0.749019622802734]);
plot3 = plot(data.time,data.y_ftx(6,:)','lineWidth',3.0,'Parent',axes1);hold on;
set(plot3,'color',[0.929411768913269 0.694117665290833 0.125490203499794]);

%plot(RNEA_t,RNEA_ftx(:,4:6),'--');
%set(h,'color',[0.5 0.5 0.5]);
leg = legend('$\tau_x$','$\tau_y$','$\tau_z$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Torque [Nm]','FontSize',20);
axis tight; 


%% plot torque estimates

fig = figure(7);
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');

%plot(data.time,res.tau_foot,'r',data.time,res.tau_leg,'b','lineWidth',2.0);
plot1 = plot(data.time,res.tau_foot,'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,res.tau_leg,'lineWidth',3.0, 'Parent',axes1);hold on;
set(plot2,'color',[0 0.498039215803146 0]);
leg = legend('$\tau_1$','$\tau_2$');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Torque [Nm]','FontSize',20);
axis tight; hold on;   


%% plot torque variance
% 
fig = figure(8);
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');

%figure(8);
% axes1 = axes('Parent',fig,'FontSize',16);
% box(axes1,'on');
% hold(axes1,'on');
% plot1 = plot(axes1,data.time,res.tau_foot, 'r', 'LineWidth', 3.0);
% %plot(axes_handle, time,avg, 'r', 'LineWidth', width);
%  hold on;
%  stddev = squeeze(2*sqrt(res.Stau_foot));
%  size(res.tau_foot)
%  minVals = res.tau_foot-stddev';
%  maxVals = res.tau_foot+stddev';
% %  ymin = min(minVals);
% %  ymax = max(maxVals);
%  timeCont = [data.time; fliplr(data.time)];
%  stdArea = [minVals', fliplr(maxVals')];
%  size(timeCont)
%  size(stdArea)
%  fHandle = fill(timeCont', stdArea, 'b');
%  alpha(fHandle, 0.5);
% 
% 
lineProps = {'LineWidth', 3.0};
shadedErrorBar(data.time,res.tau_foot,squeeze(2*sqrt(res.Stau_foot)), lineProps);hold on;
shadedErrorBar(data.time,res.tau_leg,squeeze(2*sqrt(res.Stau_foot)),lineProps);
% 
% %plot(tempT,res.tau_foot,'b',tempT,res.tau_leg,'r','lineWidth',2.0);

leg = legend('$\sigma tau_1$','$\dot tau_2$');
set(leg,'Interpreter','latex');
% set(leg,'FontSize',18);
xlabel('Time [s]');
ylabel('Torque [Nm]');
axis tight; hold on;   



%% save plots
figBaseFolder = './experiments/humanFixedBase/plots';
if(exist(figBaseFolder,'dir')==0)
    mkdir(figBaseFolder);
end

plotFigBaseName = strcat(figBaseFolder,sprintf('/PNEA_human%d_trial%d_',subjectID,trialID));
selectedFigList = [1,2,3,4,5,6,7,8];
FigName = {'jointAngles',...
            'jointVelocities',...
            'imuAccelerometer',...
            'imuGyro',...
            'forcePlateForce',...
            'forcePlateTorque',...
            'torques',...
            'TorquesVariance'};

        
for i = 1:length(selectedFigList)
%      figure(selectedFigList(i))
%      set(gca,'FontSize',12);
%      set(gcf,'Renderer','OpenGL');
%      print('-dpdf','-r300',strcat(plotFigBaseName,FigName{i}),'-opengl');
     
     figureHandle= figure(selectedFigList(i));
     ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');

    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

    print(figureHandle, '-dpdf', strcat(plotFigBaseName,FigName{i}));

end