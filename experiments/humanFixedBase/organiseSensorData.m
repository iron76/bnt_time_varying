% organiseSensorData
% Script to organise the collected sensor data of the Human-Dynamics 
% estimation experiment. The sensors are the VICON markers, IMU
% placed on the chest and force place on the bottom of the foot. 
% The quantities being computed are joint angles q1, q2, qDot1, qDot2,
% the 3D linear acceleration a_imu, 3D angular velocity omega_imu, and the 6D wrench
% measurements from the force place. The IMU and forceplate information
% is organised in the senor frames.
%
% Author: Naveen Kuppuswamy (naveen.kuppuswamy@iit.it)
% iCub Facility, Istituto Italiano di Tecnologia, 21 January 2016


%% load the synchronised dataset (including VICON and IMU data)
load('./experiments/humanFixedBase/intermediateDataFiles/synchronisedSensorData.mat');

isTest = 'false';

subjectList = 1;
trialList = 1 ; 


%% iterate through each trial computing transforms each time
for subjectID = subjectList
    fprintf('\n---------\nSubject : %d\nTrial : ',subjectID);
    for trialID = trialList
        fprintf('%d, ',trialID);
        
        currentTrial = synchronisedData(subjectID,trialID);
        
        pSelec = size(currentTrial.P_G_lhee,1);
        
        %% Extracting VICON marker trajectories
        P_G_lhee = currentTrial.P_G_lhee(1:pSelec,:);
        P_G_ltoe = currentTrial.P_G_ltoe(1:pSelec,:);
        P_G_rhee = currentTrial.P_G_rhee(1:pSelec,:);
        P_G_rtoe = currentTrial.P_G_rtoe(1:pSelec,:);
        P_G_lankle = currentTrial.P_G_lankle(1:pSelec,:);
        P_G_lhip = currentTrial.P_G_lhip(1:pSelec,:);
        P_G_lsho = currentTrial.P_G_lsho(1:pSelec,:);
        P_G_rankle = currentTrial.P_G_rankle(1:pSelec,:);
        P_G_rhip = currentTrial.P_G_rhip(1:pSelec,:);
        P_G_rsho = currentTrial.P_G_rsho(1:pSelec,:);
        P_G_tors = currentTrial.P_G_tors(1:pSelec,:);
        P_G_imuA = currentTrial.P_G_imuA(1:pSelec,:);
        P_G_imuB = currentTrial.P_G_imuB(1:pSelec,:);
        P_G_imuC = currentTrial.P_G_imuC(1:pSelec,:);
        
        %% Computing q1 and q2 angles 

        % computing P_G_1, P_G_2, P_G_3 for all the time
        P_G_1 = computeCentroidOfPoints(P_G_lankle,P_G_rankle);
        P_G_2 = computeCentroidOfPoints(P_G_lhip,P_G_rhip);
        %P_G_3 = computeCentroidOfTriangle(P_G_lsho,P_G_rsho,P_G_tors);
        P_G_3 = computeCentroidOfPoints(P_G_lsho,P_G_rsho);
        
        % JOINT ANGLE q1
        len = size(P_G_1,1);
        
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
    
        
        %% Using Savitzky-Golay filtering for computing derivatives

        % to do: tune window and polyn order.
        window = 57;
        [~, diffCoeff] = sgolay_wrapper(3, window);
        %diffCoeff is a matrix of (polynomialOrder-1) columns where:
        %- ( ,1) --> coefficient for S-Golay as smoother;
        %- ( ,2) --> coefficient for S-Golay as 1st differentiator;
        %- ( ,3) --> coefficient for S-Golay as 2nd differentiator;
        %  .   
        %  .
        %  .
        %- ( ,polynomialOrder-1) --> coefficient for S-Golay as (polynomialOrder) differentiator;
        
        halfWindow  = ((window+1)/2) -1;
        dq1_sg = zeros(len, 1);
        dq2_sg = zeros(len, 1);
        ddq1_sg = zeros(len, 1);
        ddq2_sg = zeros(len, 1);
        
        for n = (window+1)/2:len-(window+1)/2,
              % 1st differential
              dq1_sg(n) = dot(diffCoeff(:,2),q1(n - halfWindow:n + halfWindow));
              dq2_sg(n) = dot(diffCoeff(:,2),q2(n - halfWindow:n + halfWindow));
              % 2nd differential
              ddq1_sg(n) = dot(diffCoeff(:,3),q1(n - halfWindow:n + halfWindow));
              ddq2_sg(n) = dot(diffCoeff(:,3),q2(n - halfWindow:n + halfWindow));
        end
        
        %% Time step 
        % depends on sampling frequency so we extract it from time stamps
        delT = mean(diff(currentTrial.t_vicon(:,1:pSelec)));
        
        dq1_sg = dq1_sg ./ delT;
        dq2_sg = dq2_sg ./ delT;
        ddq1_sg = ddq1_sg ./ (delT)^2;
        ddq2_sg = ddq2_sg ./ (delT)^2;
        
        %% plot q, dq, ddq
        
        figure;
        subplot(311);
        plot1 = plot(currentTrial.t_vicon(:,1:pSelec),q1.*(180/pi),'lineWidth',1.0); hold on;
        set(plot1,'color',[1 0 0]);
        plot2= plot(currentTrial.t_vicon(:,1:pSelec),q2.*(180/pi),'lineWidth',1.0); hold on;
        set(plot2,'color',[0 0.498039215803146 0]);
        leg = legend('$q_1$','$q_2$','Location','northeast');
        title('Joint Quantities','FontSize',15);
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',15);
        xlabel('Time [s]','FontSize',15);
        ylabel('Joint Angle [deg]','FontSize',15);
        axis tight;
        grid on;  

        subplot(312);
        plot1 = plot(currentTrial.t_vicon(:,1:pSelec),(180/pi)*dq1_sg,'lineWidth',1.0); hold on;
        set(plot1,'color',[1 0 0]);
        plot2= plot(currentTrial.t_vicon(:,1:pSelec),(180/pi)*dq2_sg,'lineWidth',1.0); hold on;
        set(plot2,'color',[0 0.498039215803146 0]);
        leg = legend('$\dot q_{1}$','$\dot q_{2}$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',15);
        xlabel('Time [s]','FontSize',15);
        ylabel('Joint Velocity [deg/s]','FontSize',15);
        axis tight;
        grid on;
        
        subplot(313);
        plot1 = plot(currentTrial.t_vicon(:,1:pSelec),(180/pi)*ddq1_sg,'lineWidth',1.0); hold on;
        set(plot1,'color',[1 0 0]);
        plot2= plot(currentTrial.t_vicon(:,1:pSelec),(180/pi)*ddq2_sg,'lineWidth',1.0); hold on;
        set(plot2,'color',[0 0.498039215803146 0]);
        leg = legend('$\ddot q_{1}$','$\ddot q_{2}$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',15);
        xlabel('Time [s]','FontSize',15);
        ylabel('Joint Acceleration [deg/s^2]','FontSize',15);
        axis tight;
        grid on;
        
       
        %% IMU data (a_imuLin, omega_imu)
      
        % Plotting raw data coming from IMU sensor in IMU frame
        
        figure;
        subplot(211);
        plot(currentTrial.t_imu,currentTrial.a_imu_imulin'); axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Linear Acceleration [m/sec^2]','FontSize',15);
        legend('$a_x$','$a_y$','$a_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',15);
        legend('a^{IMU}_x','a^{IMU}_y','a^{IMU}_z','Location','northeast');
        title('Raw IMU data of link 2 (IMU frame)','FontSize',15);
        grid on;
        
        subplot(212);
        plot(currentTrial.t_imu,currentTrial.v_imu_imurot');  axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Angular Velocity [rad/s]','FontSize',15);
        legend('$w_x$','$w_y$','$w_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',15);
        legend('w^{IMU}_x','w^{IMU}_y','w^{IMU}_z','Location','northeast');
        grid on;
        
        a_imu = zeros(size(q1,1),3);
        omega_imu = zeros(size(q1,1),3);
        
          for i = 1:length(currentTrial.t_vicon)   
                a_imu(i,:) =  currentTrial.a_imu_imulin(i,:); 
                omega_imu(i,:) =  currentTrial.v_imu_imurot(i,:); 
          end
          
        
        %% Extracting forceplate wrench (f_fp)
        
        % Plotting raw data coming from force plate sensor in sensor frame
        f_fp = currentTrial.f_fp(1:pSelec,:);
       
        figure;
        subplot(211);
        plot(currentTrial.t_vicon,f_fp(:,4:6)); axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Force [N]','Fontsize',15);
        title('Wrench measured in force plate (Fp frame)','FontSize',15);
        legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
        grid on;
        
        subplot(212);
        plot(currentTrial.t_vicon,f_fp(:,1:3)); axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Moment [Nm]','FontSize',15);
        legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
        grid on;
 
        %% Organising into a structure
        
        processedSensorData(subjectID,trialID).q1 = q1;
        processedSensorData(subjectID,trialID).q2 = q2;
        processedSensorData(subjectID,trialID).dq1 = dq1_sg;
        processedSensorData(subjectID,trialID).dq2 = dq2_sg;
        processedSensorData(subjectID,trialID).ddq1 = ddq1_sg;
        processedSensorData(subjectID,trialID).ddq2 = ddq2_sg;
        processedSensorData(subjectID,trialID).a_imu = a_imu';
        processedSensorData(subjectID,trialID).omega_imu = omega_imu';
        processedSensorData(subjectID,trialID).t = currentTrial.t_vicon;
        processedSensorData(subjectID,trialID).imu = [a_imu omega_imu]';
        processedSensorData(subjectID,trialID).f_fp = f_fp';
    end
    fprintf('\n');
end

%% save result in the preProcessingDataFiles folder
if(strcmp(isTest,'true')~=1)
    save('./experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat','processedSensorData');
end