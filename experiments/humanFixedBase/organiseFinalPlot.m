% organiseFinalPlot
% Script that organises final plot

clear; close all; clc;

figFolder = './experiments/humanFixedBase/plots';
if(exist(figFolder,'dir')==0)
    mkdir(figFolder);
end

subjectList = 1;
trialList = 1;  

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d ',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d\n',trialID);
        
  


%% load dataset

%% plot joint angles

%% plot joint velocities

%% plot IMU accelerometer

%% plot IMU gyro

%% plot forcePlate force

%% plot forcePlate torque

%% plot torque estimates

%% plot torque variance
%
% %% save plots
% figBaseFolder = './experiments/humanFixedBase/plots';
% if(exist(figBaseFolder,'dir')==0)
%     mkdir(figBaseFolder);
% end
% 
% plotFigBaseName = strcat(figBaseFolder,sprintf('/PNEA_human%d_trial%d_',subjectID,trialID));
% selectedFigList = [1,2,3,4,5,6,7,8];
% FigName = { 'jointAngles',...
%             'jointVelocities',...
%             'imuAccelerometer',...
%             'imuGyro',...
%             'forcePlateForce',...
%             'forcePlateTorque',...
%             'torques',...
%             'TorquesVariance'};
% 
%         
% for i = 1:length(selectedFigList)
% %      figure(selectedFigList(i))
% %      set(gca,'FontSize',12);
% %      set(gcf,'Renderer','OpenGL');
% %      print('-dpdf','-r300',strcat(plotFigBaseName,FigName{i}),'-opengl');
%      
%      figureHandle= figure(selectedFigList(i));
%      ti = get(gca,'TightInset');
%     set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
% 
%     set(gca,'units','centimeters')
%     pos = get(gca,'Position');
%     ti = get(gca,'TightInset');
% 
%     set(gcf, 'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% 
%     print(figureHandle, '-dpdf', strcat(plotFigBaseName,FigName{i}));
% 
% end

  end
    fprintf('\n');
end