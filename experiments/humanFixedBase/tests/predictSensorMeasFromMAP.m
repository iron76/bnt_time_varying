% predictSensorMeasFromMAP
%
% 

clear;clc;close all;

%% testOptions

plotSensorFramePrediction = true; 

%% selected subjects and trials
subjectList = 1;
trialList = 1;  

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d\nTrial : ',subjectID);
    for trialID = trialList
         fprintf('%d, ',trialID);
         
    %% load data
    load('./experiments/humanFixedBase/intermediateDataFiles/BERDYFormattedSensorData.mat');
   
    currentTrial = BERDYFormattedSensorData(subjectID,trialID); 
    data = currentTrial.data;
    dataTime = currentTrial.data.dataTime;
    len = length(dataTime);
    
    q = currentTrial.data.q;
    dq = currentTrial.data.dq;
    ddq = currentTrial.data.ddq;
    
    % in ang-lin notation
    ys_linkFrame_imu = currentTrial.data.ys_linkFrame_imu';
    ys_link0Frame_fts = currentTrial.data.ys_link0Frame_fts';
    ys_link1Frame_fts = currentTrial.data.ys_link1Frame_fts';

%%
    sens.parts    = {'leg'         ,'torso'};       %force of the forceplate is trasmitted into the leg
    sens.labels   = {'fts'         ,'imu'  };
    sens.ndof     = {6             ,6      };

    label_to_plot = {'fts'         ,'imu'  };
    load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');
    
    currentModel = humanThreeLinkModelFromURDF(subjectID);
    
    humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
    humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 
   
    dmodel = currentModel.dmodel;                       %deterministic model 
    ymodel  = humanThreeLinkSens(dmodel, sens);  
   
    dmodel  = autoTreeStochastic(dmodel, 1e-5, 1e-4);   % probabilistic model for D equation (added Sv and Sw)
    ymodel  = humanThreeLinkSensStochastic(ymodel);     % probabilistic model for Y(q,dq) d = y (added Sy)
   
    myModel = model(dmodel);
    mySens  = sensors(ymodel);  
   
    myMAP  = MAP(myModel, mySens);
    
    %% Computing tau using Newton-Euler with Featherstone toolbox

    v_2_2 = zeros(size(q,1),6);     % spatial velocity link2
    fx = zeros (6,1);
    
    fext    = cell(1,2);
    for i = 1 : dmodel.NB
         fext{i}    = fx;
    end

    for i = 1:len
        
         [tau_i, a_i, v_i, fB_i, f_i] = IDv(dmodel, q(i,:), dq(i,:), ddq(i,:), fext);
         v_2_2(i,:) = v_i{2}';
    end

    %%  ========================= Sensor prediction ===========================

    load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');

    currentTrialSens = sensorLinkTransforms(subjectID,trialID);
    
    X_imu_2 = currentTrialSens.X_imu_2;
    XStar_fp_0 = currentTrialSens.XStar_fp_0;
    XStar_0_1 = currentTrialSens.XStar_0_1;
    
    
    load('./experiments/humanFixedBase/intermediateDataFiles/MAPresults.mat');

    currentMAP = MAPresults(subjectID,trialID);
    
    d = currentMAP.MAPres.d;
    b_Y = currentMAP.b_Y;
    %Ymatrix = currentMAP.MAPres.Ymatrix;
        
    % forcing myMAP.IDsens.sensorsParams.Y to be equal Ymatrix
    
    
  
    y_pred_MAP = zeros (myMAP.IDmeas.m,len);
    for i =1:len
         %myMAP.IDsens.sensorsParams.Y = Ymatrix{i};
         y_pred_MAP(:,i) = myMAP.simY(d(:,i));              % without b_Y
         %y_pred_MAP(:,i) = y_pred_MAP(:,i) + b_Y(:,i);     % adding b_Y
    end   
    
    %% Plot predictions 
    if(plotSensorFramePrediction)
       %% Acceleration prediction
       
        fig = figure();
        % fig = figure('name','xxx');
        axes1 = axes('Parent',fig,'FontSize',16);
        box(axes1,'on');
        hold(axes1,'on');
        grid on;
       
        subplot(221)
        plot(dataTime,y_pred_MAP(),dataTime,y_pred_MAP(),dataTime,y_pred_MAP(), 'lineWidth',2.0);
        leg = legend('$a_x$','$a_y$','$a_z$','Location','southeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Prediction ','FontSize',15);
        title('Linear acceleration [m/sec^2]','FontSize',14);
        axis tight;
        grid on; 

        subplot(223)
        plot(dataTime,imu(:,1),dataTime,imu(:,2),dataTime,imu(:,3),'lineWidth',2.0);
        leg = legend('$a_x$','$a_y$','$a_z$','Location','southeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Actual','FontSize',15);
        axis tight;
        grid on; 

        %% Gyroscope prediction comparison

        subplot(222)
        plot(dataTime,omega_imu_imuPred(1,:),dataTime,omega_imu_imuPred(2,:),dataTime,omega_imu_imuPred(3,:), 'lineWidth',2.0);
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

        
        
        
        %% Force prediction 
    
        fig = figure();
        % fig = figure('name','xxx');
        axes1 = axes('Parent',fig,'FontSize',16);
        box(axes1,'on');
        hold(axes1,'on');
        grid on;
        
        subplot(221);
        plot(dataTime,y_pred_MAP(1,:),dataTime,y_pred_MAP(2,:),dataTime,y_pred_MAP(3,:), 'lineWidth',2.0);
        leg = legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Prediction','FontSize',15);
        title('Force f_1[N] in frame associated to link1 ','FontSize',15);
        hold on; 
        axis tight;
        grid on;

        subplot(223)
        plot(dataTime, ys_link1Frame_fts(4,:),dataTime, ys_link1Frame_fts(5,:),dataTime, ys_link1Frame_fts(6,:),'lineWidth',2.0);
        leg = legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Actual','FontSize',15);
        axis tight;
        grid on;
        
        %% Moment prediction 

        subplot(222);
        plot(dataTime,y_pred_MAP(4,:),dataTime,y_pred_MAP(5,:),dataTime,y_pred_MAP(6,:), 'lineWidth',2.0);
        leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Prediction','FontSize',15);
        title('Moment_1 [Nm] in frame associated to link1','FontSize',15);
        hold on;
        axis tight;
        grid on;

        subplot(224)
        plot(dataTime, ys_link1Frame_fts(1,:),dataTime, ys_link1Frame_fts(2,:),dataTime, ys_link1Frame_fts(3,:),'lineWidth',2.0);
        leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',15);
        ylabel('Actual','FontSize',15);
        axis tight;
        grid on;
        
        
        
        
    end    


%         axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
%         text(0.5, 0.99,'\bf Prediction comparison (IMU frame)','HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',14);
% 
%         %% Force prediction comparison 
%     
%         fig = figure();
%         % fig = figure('name','xxx');
%         axes1 = axes('Parent',fig,'FontSize',16);
%         box(axes1,'on');
%         hold(axes1,'on');
%         grid on;
% 
%         subplot(221);
%         plot(dataTime,f_fp_fpPred(4,:),dataTime,f_fp_fpPred(5,:),dataTime,f_fp_fpPred(6,:), 'lineWidth',2.0);
%         leg = legend('$F_x$','$F_y$','$F_z$','Location','northeast');
%         set(leg,'Interpreter','latex');
%         set(leg,'FontSize',18);
%         xlabel('Time [s]','FontSize',15);
%         ylabel('Prediction','FontSize',15);
%         title('Force [N]','FontSize',15);
%         hold on; 
%         axis tight;
%         grid on;
% 
%         subplot(223)
%         plot(dataTime,wrench_fp_fp(:,4),dataTime,wrench_fp_fp(:,5),dataTime,wrench_fp_fp(:,6),'lineWidth',2.0);
%         leg = legend('$F_x$','$F_y$','$F_z$','Location','northeast');
%         set(leg,'Interpreter','latex');
%         set(leg,'FontSize',18);
%         xlabel('Time [s]','FontSize',15);
%         ylabel('Actual','FontSize',15);
%         axis tight;
%         grid on;
% 
%         %% Moment prediction comparison 
% 
%         subplot(222);
%         plot(dataTime,f_fp_fpPred(1,:),dataTime,f_fp_fpPred(2,:),dataTime,f_fp_fpPred(3,:), 'lineWidth',2.0);
%         leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
%         set(leg,'Interpreter','latex');
%         set(leg,'FontSize',18);
%         xlabel('Time [s]','FontSize',15);
%         ylabel('Prediction','FontSize',15);
%         title('Moment [Nm]','FontSize',15);
%         hold on;
%         axis tight;
%         grid on;
% 
%         subplot(224)
%         plot(dataTime,wrench_fp_fp(:,1),dataTime,wrench_fp_fp(:,2),dataTime,wrench_fp_fp(:,3),'lineWidth',2.0);
%         leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
%         set(leg,'Interpreter','latex');
%         set(leg,'FontSize',18);
%         xlabel('Time [s]','FontSize',15);
%         ylabel('Actual','FontSize',15);
%         axis tight;
%         grid on;
%     
%         axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
%         text(0.5, 0.99,'\bf Prediction comparison (Fp frame)','HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',14);
%         
%         
   
    end
      fprintf('\n');
end

