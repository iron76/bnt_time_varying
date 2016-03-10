% plotSynchronisedSensorData

clc; clear; close all;

%% testOptions
plotIMUData = true;
plotForcePlateData = true;

%% load the synchronised dataset
load('./experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat');

isTest = 'false';

subjectList = 1;
trialList = 1 ; 

%% iterate through each trial 
for subjectID = subjectList
    fprintf('\n---------\nSubject : %d ',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d ',trialID);
        
        currentTrial = synchronisedData(subjectID,trialID);
        dataTime = currentTrial.dataTime;
        
        aLin_imu_imu = currentTrial.aLin_imu_imu;
        omega_imu_imu = currentTrial.omega_imu_imu;
        wrench_fp_fp = currentTrial.wrench_fp_fp;

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
        
    end
    fprintf('\n');
end
