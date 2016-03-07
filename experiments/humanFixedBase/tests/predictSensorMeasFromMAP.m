% predictSensorMeasFromMAP
% Script makes a prediction of the sensor measurements simulating output of
% MAP class. Sensors are the VICON markers, IMU placed on the chest 
% and force place on the bottom of the foot.
%
% Toolbox requirements:
% - iDynTree - mex
% - Featherstone toolbox (v2)

clear;clc;close all;

%% testOptions

plotSensorPrediction = true; 
plotResultsVariances = true;

%% selected subjects and trials
subjectList = 1;
trialList = 1;  

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d',subjectID);
    for trialID = trialList
         fprintf('\nTrial : %d\n ',trialID);
         
    %% Load data 

    load('./experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat');

    currentTrial = processedSensorData(subjectID,trialID);
        
    q1 = currentTrial.q1;
    q2 = currentTrial.q2; 
    dq1 = currentTrial.dq1;
    dq2 = currentTrial.dq2;
    ddq1 = currentTrial.ddq1;
    ddq2 = currentTrial.ddq2;    
    dataTime = currentTrial.dataTime;
    
    len = length(dataTime);
    q = [q1 q2];
    dq = [dq1 dq2];
    ddq = [ddq1 ddq2];
        
    wrench_fp_fp = currentTrial.wrench_fp_fp;
    imu = currentTrial.imu;

    %% Load models and build structure
    
    sens.parts    = {'leg'         ,'torso'};           %force of the forceplate is trasmitted into the leg
    sens.labels   = {'fts'         ,'imu'  };
    sens.ndof     = {6             ,6      };

    label_to_plot = {'fts'         ,'imu'  };
    load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');
    
    currentModel = humanThreeLinkModelFromURDF(subjectID);
    
    humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
    humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 
   
    dmodel = currentModel.dmodel;                       % deterministic model 
    ymodel  = humanThreeLinkSens(dmodel, sens);  
   
    dmodel  = autoTreeStochastic(dmodel, 1e-6);   % probabilistic model for D equation (added Sv and Sw)
    ymodel  = humanThreeLinkSensStochastic(ymodel);     % probabilistic model for Y(q,dq) d = y (added Sy)
   
    myModel = model(dmodel);
    mySens  = sensors(ymodel);  
   
    myMAP  = MAP(myModel, mySens);
    
    %%  ========================= Sensor prediction (NO MAP) ===========================

%     load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');
% 
%     currentTrialSens = sensorLinkTransforms(subjectID,trialID);
%     
%     X_imu_2 = currentTrialSens.X_imu_2;
%     XStar_fp_0 = currentTrialSens.XStar_fp_0;
%     XStar_0_1 = currentTrialSens.XStar_0_1;
%     
%     
%     load('./experiments/humanFixedBase/intermediateDataFiles/MAPresults.mat');
%  
%     currentMAP = MAPresults(subjectID,trialID);
%     
%     d = currentMAP.MAPres.d;
%     b_Y = currentMAP.MAPres.b_Y;
%     Ymatrix = currentMAP.MAPres.Ymatrix;
%     
%     fprintf('Predicting sensor measurement using MAP computation\n');
% 
%     y_pred_MAP = zeros (myMAP.IDmeas.m,len);
%     for i =1:len
%          y_pred_MAP(:,i) = myMAP.simY(d(:,i));             % without b_Y
%          %y_pred_MAP(:,i) = y_pred_MAP(:,i) + b_Y(:,i);     % adding b_Y
%     end   
    
    %%  ========================= Sensor prediction ===========================
    
    load('./experiments/humanFixedBase/intermediateDataFiles/MAPresults.mat');
 
    currentMAP = MAPresults(subjectID,trialID);
    
    b_Y = currentMAP.MAPres.b_Y;
    Ymatrix = currentMAP.MAPres.Ymatrix;
    data.y = currentMAP.data.y;
    data.Sy = currentMAP.data.Sy;
    
    y_pred = zeros (myMAP.IDmeas.m,len);
    
    % Computing MAP method
    fprintf('Predicting sensor measurement using MAP computation\n');

    for i = 1 : len

        myMAP = myMAP.setState(q(i,:)',q(i,:)');
        myMAP = myMAP.setY(data.y(:,i));
        myMAP = myMAP.setYmatrix(Ymatrix{i});
        myMAP = myMAP.setBias(b_Y(:,i));
        myMAP = myMAP.solveID();
  
        res.d(:,i)       = myMAP.d;
        res.Sd(:,:,i)    = myMAP.Sd; %full() passing from sparse to double matrix
        res.Ymatrix{i,1} = myMAP.IDsens.sensorsParams.Y; 
        res.b_Y(:,i)     = myMAP.IDsens.sensorsParams.bias;
        
        %simulate output usinf simY
        res.y_pred(:,i)  = myMAP.simY(res.d(:,i));
    
         if mod(i-1,100) == 0
                fprintf('Processing %d %% of the dataset\n', round(i/len*100));
         end
    end
    
% ========end MAP 

    %% Plot predictions 
    if(plotSensorPrediction)

       %% Accelerometer prediction comparison
       
        %fig = figure();
        fig = figure('name','MAP');
        axes1 = axes('Parent',fig,'FontSize',16);
        box(axes1,'on');
        hold(axes1,'on');
        grid on;
       
        subplot(221)
        plot(dataTime,res.y_pred(10,:),dataTime,res.y_pred(11,:),dataTime,res.y_pred(12,:),'lineWidth',2.0);
        leg = legend('$a_x$','$a_y$','$a_z$','Location','southeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Prediction ','FontSize',15);
        title('Linear acceleration [m/sec^2]','FontSize',14);
        axis tight;
        grid on; 

        subplot(223)
        plot2 = plot(dataTime,imu(:,1),dataTime,imu(:,2),dataTime,imu(:,3),'lineWidth',2.0);
        leg = legend('$a_x$','$a_y$','$a_z$','Location','southeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Actual','FontSize',15);
        axis tight;
        grid on; 

        %% Gyroscope prediction comparison

        subplot(222)
        plot(dataTime,res.y_pred(7,:),dataTime,res.y_pred(8,:),dataTime,res.y_pred(9,:), 'lineWidth',2.0);
        leg = legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','southeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Prediction','FontSize',15);
        title('Angular velocity [rad/sec]','FontSize',14);
        axis tight;
        grid on; 

        subplot(224)
        plot(dataTime,imu(:,4),dataTime,imu(:,5),dataTime,imu(:,6),'lineWidth',2.0);
        leg = legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','southeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Actual','FontSize',15);
        axis tight;
        grid on; 

        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
        text(0.5, 0.99,(sprintf('MAP measurements prediction, sensor frame (Subject: %d, Trial: %d)',subjectID, trialID)),'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',14);

        %% Force prediction comparison  
    
        %fig = figure();
        fig = figure('name','MAP');
        axes1 = axes('Parent',fig,'FontSize',16);
        box(axes1,'on');
        hold(axes1,'on');
        grid on;
        
        subplot(221);
        plot(dataTime,res.y_pred(4,:),dataTime,res.y_pred(5,:),dataTime,res.y_pred(6,:), 'lineWidth',2.0);
        leg = legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Prediction','FontSize',15);
        title('Force [N]','FontSize',15);
        hold on; 
        axis tight;
        grid on;

        subplot(223)
        plot(dataTime,wrench_fp_fp(:,4),dataTime,wrench_fp_fp(:,5),dataTime,wrench_fp_fp(:,6),'lineWidth',2.0);
        leg = legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Actual','FontSize',15);
        axis tight;
        grid on;
        
        %% Moment prediction comparison

        subplot(222);
        plot(dataTime,res.y_pred(1,:),dataTime,res.y_pred(2,:),dataTime,res.y_pred(3,:), 'lineWidth',2.0);
        leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Prediction','FontSize',15);
        title('Moment [Nm]','FontSize',15);
        hold on;
        axis tight;
        grid on;

        subplot(224)
        plot(dataTime,wrench_fp_fp(:,1),dataTime,wrench_fp_fp(:,2),dataTime,wrench_fp_fp(:,3),'lineWidth',2.0);
        leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Actual','FontSize',15);
        axis tight;
        grid on;
        
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
        text(0.5, 0.99,(sprintf('MAP measurements prediction, sensor frame (Subject: %d, Trial: %d)',subjectID, trialID)),'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',14); 
        
    end    
    
    %% Plot variances
    if(plotResultsVariances)
         
        for i = 1:myMAP.IDsens.sensorsParams.m
        
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;
      
            lineProps = {'LineWidth', 3.0};
            shad1 = shadedErrorBar(dataTime,data.y(i,:),sqrt(data.Sy(i,:)),lineProps);
            plot1= plot(dataTime,res.y_pred(i,:),'lineWidth',1.0);hold on;
            set(plot1,'color',[1 0 0]);
            xlabel('Time [s]','FontSize',15);
            leg = legend([shad1.mainLine,shad1.patch, plot1], {'sensorData','dataVariance','predictionData'},'Location','southeast');
            title('Comparison between real data and MAP prediction data','FontSize',13);

            if (i >=1 && i<4) 
                ylabel('f1: Moment [Nm]','FontSize',15);
            elseif (i >=4 && i<7) 
                ylabel('f1: Force[N]','FontSize',15);
            elseif (i >=7 && i<10)
                ylabel('a2: Angular acceleration [rad/s]','FontSize',15);
            elseif (i >=10 && i<13)
                ylabel('a2: Linear acceleration [m/s^2]','FontSize',15);
            elseif (i >=13 && i<19)
                ylabel('External force ftx1 [N]','FontSize',15);
            elseif (i >=19 && i<25)
                ylabel('External force ftx2 [N]','FontSize',15);
            elseif (i ==25)
                ylabel('Joint acceleration ddq1 [m/s^2]','FontSize',15);
            elseif (i ==26)
                ylabel('Joint acceleration ddq2 [m/s^2]','FontSize',15);
            end
            
            axis tight;
            grid on; 
        end    
    end 
   end
      fprintf('\n');
end
