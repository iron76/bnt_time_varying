close all;

% subjectID = 1;
% trialID = 1;

%% load dataset
%load(sprintf('./experiments/humanFixedBase/data/savedBERDYresult_subj%d_trial%d.mat',subjectID,trialID),'res','data','myPNEA');

%% plot joint angles


%% plot joint velocities


%% plot IMU accelerometer
%[RNEA_imu,RNEA_ftx,RNEA_t] = sensorDataConsistencyCheck(subjectID,trialID,'noplots');


%% plot IMU gyro


%% plot forcePlate force


%% plot forcePlate torque


%% plot torque estimates



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