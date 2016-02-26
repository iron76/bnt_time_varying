% organiseSensorData
% Script to organise the collected sensor data.  
% Computed clustered quantities are:
% - joint angles q1 and q2;
% - joint velocities dq1 and dq2;
% - joint accelerations ddq1 and ddq2;
% Measured clustered quantities are:
% - linear acceleration from IMU a_imu;
% - angular velocity from IMU omega_imu;
% - wrench from force plate f_fp_fp;
% Measured quantities are in sensor frame.

clc; clear; close all;

%% testOptions
plotJointQuantities = false;
plotIMUData = false;
plotForcePlateData = false;

%% load the synchronised dataset
load('./experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat');

isTest = 'false';

subjectList = 1:12;
trialList = 1:4 ; 


%% iterate through each trial 
for subjectID = subjectList
    fprintf('\n---------\nSubject : %d ',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d ',trialID);
        
        currentTrial = synchronisedData(subjectID,trialID);
        
        P_G_lhee = currentTrial.P_G_lhee;
        P_G_ltoe = currentTrial.P_G_ltoe;
        P_G_rhee = currentTrial.P_G_rhee;
        P_G_rtoe = currentTrial.P_G_rtoe;
        P_G_lankle = currentTrial.P_G_lankle;
        P_G_lhip = currentTrial.P_G_lhip;
        P_G_lsho = currentTrial.P_G_lsho;
        P_G_rankle = currentTrial.P_G_rankle;
        P_G_rhip = currentTrial.P_G_rhip;
        P_G_rsho = currentTrial.P_G_rsho;
        P_G_tors = currentTrial.P_G_tors;
        P_G_imuA = currentTrial.P_G_imuA;
        P_G_imuB = currentTrial.P_G_imuB;
        P_G_imuC = currentTrial.P_G_imuC;
        
        samplingTime = currentTrial.samplingTime;
        dataTime = currentTrial.dataTime;
        len = length(P_G_lhee);
        
        aLin_imu_imu = currentTrial.aLin_imu_imu;
        omega_imu_imu = currentTrial.omega_imu_imu;
        wrench_fp_fp = currentTrial.wrench_fp_fp;
        
        %% computing q1 and q2 angles 

        P_G_1 = computeCentroidOfPoints(P_G_lankle,P_G_rankle);
        P_G_2 = computeCentroidOfPoints(P_G_lhip,P_G_rhip);
        P_G_3 = computeCentroidOfPoints(P_G_lsho,P_G_rsho);
        
        % JOINT ANGLE q1
        l1 = (P_G_2 - P_G_1);
        q1 = zeros (len, 1);
        
        for i = 1 : len;
            q1(i) =atan2(-l1(i,2),l1(i,3));
        end
       
        % JOINT ANGLE q2
        l2 = (P_G_3-P_G_2);
        q_temp = zeros (len, 1);

        for i = 1 : len;
            q_temp(i) = atan2(-l2(i,2),l2(i,3));
        end
        
        q2 = q_temp-q1;
          
        %% computing dq1,dq2,ddq1,ddq2

        [dq1,ddq1] = SgolayDerivation(3,57,q1,samplingTime);
        [dq2,ddq2] = SgolayDerivation(3,57,q2,samplingTime);
        
        %% plot q, dq, ddq
        
        if(plotJointQuantities)
        
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;
        
            subplot(311);
            plot1 = plot(dataTime,q1.*(180/pi),'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(dataTime,q2.*(180/pi),'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$q_1$','$q_2$','Location','northeast');
            %title('Joint Quantities','FontSize',15);
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Angle [deg]','FontSize',15);
            axis tight;
            grid on;  
            title(sprintf('Subject %d, Trial %d, Joint Quantities',subjectID,trialID));
        
            subplot(312);
            plot1 = plot(dataTime,(180/pi)*dq1,'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(dataTime,(180/pi)*dq2,'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$\dot q_{1}$','$\dot q_{2}$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Velocity [deg/s]','FontSize',15);
            axis tight;
            grid on;
        
            subplot(313);
            plot1 = plot(dataTime,(180/pi)*ddq1,'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(dataTime,(180/pi)*ddq2,'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$\ddot q_{1}$','$\ddot q_{2}$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Acceleration [deg/s^2]','FontSize',15);
            axis tight;
            grid on;
        end
        
        %% IMU data (a_imuLin, omega_imu)
      
        if(plotIMUData)
            % Plotting raw data coming from IMU sensor in IMU frame     
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on; 
        
            subplot(211);
            plot(dataTime, aLin_imu_imu); axis tight;
            leg = legend('$a_x$','$a_y$','$a_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Linear Acceleration [m/sec^2]','FontSize',15);
            %title('Raw IMU data of link 2 (IMU frame)','FontSize',15);
            grid on;
            title(sprintf('Subject %d, Trial %d, Raw IMU data of link 2 (IMU frame)',subjectID,trialID));
            
            subplot(212);
            plot(dataTime,omega_imu_imu);  axis tight;
            leg = legend('$w_x$','$w_y$','$w_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Angular Velocity [rad/s]','FontSize',15);
            grid on;
        end
  
        %% forceplate wrench data (f_fp_fp)
        
        if(plotForcePlateData)
            % Plotting raw data coming from force plate sensor in sensor frame      
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;
        
            subplot(211);
            plot(dataTime,wrench_fp_fp(:,4:6)); axis tight;
            leg =legend('$F_x$','$F_y$','$F_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Force [N]','Fontsize',15);
            %title('Wrench measured in force plate (Fp frame)','FontSize',15);
            grid on;
            title(sprintf('Subject %d, Trial %d, Wrench measured in force plate (Fp frame)',subjectID,trialID));
         
            subplot(212);
            plot(dataTime,wrench_fp_fp(:,1:3)); axis tight;
            leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Moment [Nm]','FontSize',15);
            grid on;
        end
 
        %% data storing
        
        processedSensorData(subjectID,trialID).samplingTime = samplingTime;
        processedSensorData(subjectID,trialID).q1 = q1;
        processedSensorData(subjectID,trialID).q2 = q2;
        processedSensorData(subjectID,trialID).dq1 = dq1;
        processedSensorData(subjectID,trialID).dq2 = dq2;
        processedSensorData(subjectID,trialID).ddq1 = ddq1;
        processedSensorData(subjectID,trialID).ddq2 = ddq2;
        processedSensorData(subjectID,trialID).dataTime = dataTime;
        processedSensorData(subjectID,trialID).imu = [aLin_imu_imu omega_imu_imu];
        processedSensorData(subjectID,trialID).wrench_fp_fp = wrench_fp_fp;
        
        processedSensorData(subjectID,trialID).P_G_ltoe =P_G_ltoe;
        processedSensorData(subjectID,trialID).P_G_lhee =P_G_lhee;
        processedSensorData(subjectID,trialID).P_G_lankle =P_G_lankle;
        processedSensorData(subjectID,trialID).P_G_lhip =P_G_lhip;
        processedSensorData(subjectID,trialID).P_G_lsho =P_G_lsho; 
        
        processedSensorData(subjectID,trialID).P_G_rtoe =P_G_rtoe;
        processedSensorData(subjectID,trialID).P_G_rhee = P_G_rhee;
        processedSensorData(subjectID,trialID).P_G_rankle =P_G_rankle;
        processedSensorData(subjectID,trialID).P_G_rhip = P_G_rhip;
        processedSensorData(subjectID,trialID).P_G_rsho = P_G_rsho;
        
        processedSensorData(subjectID,trialID).P_G_tors =P_G_tors;
        processedSensorData(subjectID,trialID).P_G_imuA = P_G_imuA;
        processedSensorData(subjectID,trialID).P_G_imuB =P_G_imuB;
        processedSensorData(subjectID,trialID).P_G_imuC =P_G_imuC;
        
        processedSensorData(subjectID,trialID).P_G_1 = P_G_1;
        processedSensorData(subjectID,trialID).P_G_2 = P_G_2;
        processedSensorData(subjectID,trialID).P_G_3 = P_G_3;
        
    end
    fprintf('\n');
end

%% save result in the preProcessingDataFiles folder
if(strcmp(isTest,'true')~=1)
    save('./experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat','processedSensorData');
end