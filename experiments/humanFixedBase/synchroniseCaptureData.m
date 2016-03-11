% synchroniseCaptureData
% Script to synchronise the captured IMU and VICON/force plate data. The 
% It looks for the spikes in beginning and end of both the sources 
% (corresponding to the little tiptoe gesture performed by subjects in the 
% start and conclusion of the bowing motion experiment). It:
%
% 1. converts eventually data in ISU;
% 2. interpolates data to obtain the same samplingTime for all of them;
% 3. computing q1,q2,dq1,dq2,ddq1,ddq2 using data coming from markers;
% 4. syncrhonises data and joint quantities previously computed:
%    4.1 re-aligns samples on the first peak cutting data before it;
%    4.2 cuts samples after the second peak;
%    4.3 cuts for both peaks  x seconds of acquisition; 
% 5. filters forceplate and IMU data;
% 6. stores data into a file.mat


% Assumptions:
% - C3D processing is used to generate the Vicon/Force plate data in .mat
% files;
% - previous processing to generate the IMU data in .mat file;


clc; clear; close all;

%% testOptions
plotInterpolatedFilteredData = true;
plotFilteredData = true;

plotFirstCutData = true;
plotJointQuantitiesFirstCut = true;

plotSecondCutData = true;
plotJointQuantitiesSecondCut = true;

plotSettlingTimeCutData = true;
plotJointQuantitiesSettlingTimeCut = true;

%% load data sources
load('./experiments/humanFixedBase/data/VICONsaveDataGen16.mat');
load('./experiments/humanFixedBase/data/imuExtractedDataGen16.mat');
addpath('./experiments/humanFixedBase/helperFunctions/');

subjectList = 7;
trialList = 1;

% Setting parameters
samplingTime = 1e-2;  % uniform sampling rate in the interpolation
settlingCutTime = 0.8;  % secs  -->TO BE SET PROPERLY! MAKE IT PARAMETRIC!
settlingCutTimeIndex = ceil(settlingCutTime/samplingTime); %ensuring integer variable

% from sensor acquisition
forcePlateSamplingTime = 1e-3;
imuSamplingTime = 1e-2;
markersSamplingTime = 1e-2;

% for filtering in frequency
cutOffFreqForcePlate = 20; %20Hz
cutOffFreqIMU = 20;  %20Hz


%% iterate through each trial 
        
 for subjectID = subjectList
    fprintf('\n---------\nSubject : %d ',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d ',trialID);
        
        %% 1. data conversion
        
        % ====VICON MARKERS
        % acquired data in mm --> converted in m
        convConst = 1e-3;
        
        P_G_ltoe_raw = convConst * subjectData(subjectID,trialID).markers.ltoe;
        P_G_lhee_raw = convConst * subjectData(subjectID,trialID).markers.lhee;
        P_G_lankle_raw = convConst * subjectData(subjectID,trialID).markers.lankle;
        P_G_lhip_raw = convConst * subjectData(subjectID,trialID).markers.lhip;
        P_G_lsho_raw = convConst * subjectData(subjectID,trialID).markers.lsho; 
        
        P_G_rtoe_raw =  convConst * subjectData(subjectID,trialID).markers.rtoe;
        P_G_rhee_raw =  convConst * subjectData(subjectID,trialID).markers.rhee;
        P_G_rankle_raw = convConst * subjectData(subjectID,trialID).markers.rankle;
        P_G_rhip_raw =  convConst * subjectData(subjectID,trialID).markers.rhip;
        P_G_rsho_raw =  convConst * subjectData(subjectID,trialID).markers.rsho; 
        
        P_G_tors_raw =  convConst * subjectData(subjectID,trialID).markers.tors; 
        P_G_imuA_raw =  convConst * subjectData(subjectID,trialID).markers.imuA;
        P_G_imuB_raw =  convConst * subjectData(subjectID,trialID).markers.imuB;
        P_G_imuC_raw =  convConst * subjectData(subjectID,trialID).markers.imuC;
        
        % ====FORCE PLATE
        % acquired force data in N --> no conversion;
        % acquired moment data in Nmm --> converted in Nm;
       
        f_raw = subjectData(subjectID,trialID).analogsFOR;
        mom_raw = convConst *  subjectData(subjectID,trialID).analogsMOM;
        
        % ====IMU
        % acquired acc data in m/s^2 --> no conversion;
        % acquired omega data in rad/s^2 --> no conversion;
        accl_raw = imuData(subjectID,trialID).accln;
        omega_raw = imuData(subjectID,trialID).gyro;
        
        %% 2. data interpolation
        
        % ====FORCE PLATE 
        t_f_raw = 0:forcePlateSamplingTime:(forcePlateSamplingTime*(size(f_raw,1) -1));
        t_f = 0:samplingTime:t_f_raw(end);
        f = interp1(t_f_raw,f_raw,t_f);
        mom = interp1(t_f_raw,mom_raw,t_f);
        
        % ====IMU      
        if(~isempty(find(diff(imuData(subjectID,trialID).t)<=0)))
            ttemp = imuData(subjectID,trialID).t;
            t_imu_raw = linspace(0,ttemp(end),length(ttemp))./1000;
        else
            t_imu_raw = imuData(subjectID,trialID).t./1000;
        end
        
        t_imu = 0:samplingTime:t_imu_raw(end);
        accl = interp1(t_imu_raw,accl_raw,t_imu);
        omega = interp1(t_imu_raw,omega_raw,t_imu);
        
        % ====VICON MARKERS
        t_markers_raw = 0:markersSamplingTime:(markersSamplingTime*size(subjectData(subjectID,trialID).markers.ltoe,1));
        t_markers_raw = t_markers_raw(1:end-1);
        
        P_G_ltoe = interp1( t_markers_raw,P_G_ltoe_raw,t_f);
        P_G_lhee = interp1( t_markers_raw,P_G_lhee_raw,t_f);
        P_G_lankle = interp1( t_markers_raw,P_G_lankle_raw,t_f);
        P_G_lhip = interp1( t_markers_raw,P_G_lhip_raw,t_f);
        P_G_lsho =interp1( t_markers_raw,P_G_lsho_raw,t_f); 
        
        P_G_rtoe = interp1(t_markers_raw,P_G_rtoe_raw,t_f);
        P_G_rhee = interp1( t_markers_raw,P_G_rhee_raw,t_f);
        P_G_rankle =interp1(t_markers_raw,P_G_rankle_raw,t_f);
        P_G_rhip = interp1( t_markers_raw,P_G_rhip_raw,t_f);
        P_G_rsho = interp1(t_markers_raw,P_G_rsho_raw,t_f); 
        
        P_G_tors =interp1( t_markers_raw,P_G_tors_raw,t_f); 
        P_G_imuA = interp1( t_markers_raw,P_G_imuA_raw,t_f);
        P_G_imuB = interp1( t_markers_raw,P_G_imuB_raw,t_f);
        P_G_imuC =interp1( t_markers_raw,P_G_imuC_raw,t_f);
        
        % filter only markers data as they are used to compute q,dq,ddq
        
        P_G_ltoe_filt = frequencyFilterSignal (P_G_ltoe, cutOffFreqForcePlate, 1/samplingTime);
        P_G_lhee_filt = frequencyFilterSignal (P_G_lhee, cutOffFreqForcePlate, 1/samplingTime);
        P_G_lankle_filt = frequencyFilterSignal (P_G_lankle, cutOffFreqForcePlate, 1/samplingTime);
        P_G_lhip_filt = frequencyFilterSignal (P_G_lhip, cutOffFreqForcePlate, 1/samplingTime);
        P_G_lsho_filt = frequencyFilterSignal (P_G_lsho, cutOffFreqForcePlate, 1/samplingTime);
        
        P_G_rtoe_filt = frequencyFilterSignal (P_G_rtoe, cutOffFreqForcePlate, 1/samplingTime);
        P_G_rhee_filt = frequencyFilterSignal (P_G_rhee, cutOffFreqForcePlate, 1/samplingTime);
        P_G_rankle_filt = frequencyFilterSignal (P_G_rankle, cutOffFreqForcePlate, 1/samplingTime);
        P_G_rhip_filt = frequencyFilterSignal (P_G_rhip, cutOffFreqForcePlate, 1/samplingTime);
        P_G_rsho_filt = frequencyFilterSignal (P_G_rsho, cutOffFreqForcePlate, 1/samplingTime);
        
        P_G_tors_filt = frequencyFilterSignal (P_G_tors, cutOffFreqForcePlate, 1/samplingTime);
        P_G_imuA_filt = frequencyFilterSignal (P_G_imuA, cutOffFreqForcePlate, 1/samplingTime);
        P_G_imuB_filt = frequencyFilterSignal (P_G_imuB, cutOffFreqForcePlate, 1/samplingTime);
        P_G_imuC_filt = frequencyFilterSignal (P_G_imuC, cutOffFreqForcePlate, 1/samplingTime);
 
        
        if(plotInterpolatedFilteredData)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;
            
            subplot(411);
            plot(t_f,f, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Force [N]','FontSize',12);
            axis tight;
            grid on;
            title(sprintf('Subject : %d, Trial : %d, interpolated data',subjectID,trialID));
            
            subplot(412);
            plot(t_f,mom, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Moment [Nmm]','FontSize',12);
            axis tight;
            grid on;
            
            subplot(413);
            plot(t_imu,accl, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('LinAcc [m/s^2]','FontSize',12);
            axis tight;
            grid on;
            
            subplot(414);
            plot(t_imu,omega, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('AngVel [rad/s]','FontSize',12);
            axis tight;
            grid on;
            
        end
        
        %% 3. computing q1,q2,dq1,dq2,ddq1,ddq2
        
        len = length(P_G_lhee);

        P_G_1 = computeCentroidOfPoints(P_G_lankle,P_G_rankle);
        P_G_2 = computeCentroidOfPoints(P_G_lhip,P_G_rhip);
        P_G_3 = computeCentroidOfPoints(P_G_lsho,P_G_rsho);
        
        % ====joint angle q1
        l1 = (P_G_2 - P_G_1);
        q1 = zeros (len, 1);
        
        for i = 1 : len;
            q1(i) =atan2(-l1(i,2),l1(i,3));
        end
       
        % ====joint angle q2
        l2 = (P_G_3-P_G_2);
        q_temp = zeros (len, 1);

        for i = 1 : len;
            q_temp(i) = atan2(-l2(i,2),l2(i,3));
        end
        
        q2 = q_temp-q1;

        [dq1,ddq1] = SgolayDerivation(3,57,q1,samplingTime);
        [dq2,ddq2] = SgolayDerivation(3,57,q2,samplingTime);
        
        %% 4. synchronising data
        
        % 4.1 synchronising peaks of force with imu 
        
        chosenF_ID = 3; % because force is on z axis
        [~,timeIndexToPeak1_f] = max(abs(f(1:round(end*0.5),chosenF_ID))); % max absolute 
        timeToPeak1_f = t_f(timeIndexToPeak1_f);

        [~,chosenA_ID] = max(abs(accl(1,:)));
        [val,timeIndexToPeak1_accl] = max(abs(accl(1:round(end*0.5),chosenA_ID)));
        timeToPeak1_accl = t_imu(timeIndexToPeak1_accl);


        % t=0 in f dal campione timeToPeak1_f
        fCut = f(timeIndexToPeak1_f:end,:);
        t_cut_vicon = t_f(timeIndexToPeak1_f:end) - t_f(timeIndexToPeak1_f);
        momCut = mom(timeIndexToPeak1_f:end,:);
        
        % t=o in acc dal campione timeToPeak1_accl
        accCut = accl(timeIndexToPeak1_accl:end,:);
        t_cut_imu = t_imu(timeIndexToPeak1_accl:end) - t_imu(timeIndexToPeak1_accl);
        omegaCut = omega(timeIndexToPeak1_accl:end,:);
        
        % synchro also q,dq,ddq
        q1 = q1(timeIndexToPeak1_f:end,:);
        q2 = q2(timeIndexToPeak1_f:end,:);
        dq1 = dq1(timeIndexToPeak1_f:end,:);
        dq2 = dq2(timeIndexToPeak1_f:end,:);
        ddq1 = ddq1(timeIndexToPeak1_f:end,:);
        ddq2 = ddq2(timeIndexToPeak1_f:end,:);
        
        if(plotFirstCutData)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            
            grid on;
            subplot(411);
            plot(t_cut_vicon,fCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Force [N] ','FontSize',12);
            axis tight;
            grid on; 
            title(sprintf('Subject : %d, Trial : %d, first cut',subjectID,trialID));
            
            subplot(412);
            plot(t_cut_vicon,momCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Moment [Nmm] ','FontSize',12);
            axis tight;
            grid on;
            
            subplot(413);
            plot(t_cut_imu,accCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('LinAcc [m/sec^2]','FontSize',12);
            axis tight;
            grid on;
            subplot(414);
            plot(t_cut_imu,omegaCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('AngVel [rad/s]','FontSize',12);
            axis tight;
            grid on;
        end
        
        if(plotJointQuantitiesFirstCut)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;
        
            subplot(311);
            plot1 = plot(t_cut_vicon,q1.*(180/pi),'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(t_cut_vicon,q2.*(180/pi),'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$q_1$','$q_2$','Location','northeast');
            %title('Joint Quantities','FontSize',15);
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Angle [deg]','FontSize',15);
            axis tight;
            grid on;  
            title(sprintf('Subject %d, Trial %d, Joint Quantities firstCut',subjectID,trialID));
        
            subplot(312);
            plot1 = plot(t_cut_vicon,(180/pi)*dq1,'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(t_cut_vicon,(180/pi)*dq2,'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$\dot q_{1}$','$\dot q_{2}$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Velocity [deg/s]','FontSize',15);
            axis tight;
            grid on;
        
            subplot(313);
            plot1 = plot(t_cut_vicon,(180/pi)*ddq1,'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(t_cut_vicon,(180/pi)*ddq2,'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$\ddot q_{1}$','$\ddot q_{2}$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Acceleration [deg/s^2]','FontSize',15);
            axis tight;
            grid on;
        end

        P_G_ltoe_cut    = P_G_ltoe   (timeIndexToPeak1_f:end,:);
        P_G_lhee_cut    = P_G_lhee  (timeIndexToPeak1_f:end,:);
        P_G_lankle_cut  = P_G_lankle (timeIndexToPeak1_f:end,:);
        P_G_lhip_cut    = P_G_lhip   (timeIndexToPeak1_f:end,:);
        P_G_lsho_cut    = P_G_lsho  (timeIndexToPeak1_f:end,:);
                   
        P_G_rtoe_cut    = P_G_rtoe  (timeIndexToPeak1_f:end,:);
        P_G_rhee_cut    = P_G_rhee (timeIndexToPeak1_f:end,:);
        P_G_rankle_cut  = P_G_rankle  (timeIndexToPeak1_f:end,:);
        P_G_rhip_cut    = P_G_rhip  (timeIndexToPeak1_f:end,:);
        P_G_rsho_cut    = P_G_rsho  (timeIndexToPeak1_f:end,:);
                    
        P_G_tors_cut    = P_G_tors  (timeIndexToPeak1_f:end,:);
        P_G_imuA_cut    = P_G_imuA  (timeIndexToPeak1_f:end,:);
        P_G_imuB_cut    = P_G_imuB  (timeIndexToPeak1_f:end,:);
        P_G_imuC_cut    = P_G_imuC  (timeIndexToPeak1_f:end,:);
        
        
        
        % 4.2 cutting data after last peak 
        
        halfIndex = round(length(fCut)*0.5);
        [~,timeIndexToPeak2_f] = max(abs(fCut(halfIndex:end,chosenF_ID))); 
        timeIndexToPeak2_f = halfIndex + timeIndexToPeak2_f;
        
        fCut = fCut(1:timeIndexToPeak2_f,:);
        t_cut_vicon = t_cut_vicon(1:timeIndexToPeak2_f);
        momCut = momCut(1:timeIndexToPeak2_f,:);
        
        P_G_ltoe_cut    = P_G_ltoe_cut   (1:timeIndexToPeak2_f,:);
        P_G_lhee_cut    = P_G_lhee_cut   (1:timeIndexToPeak2_f,:);
        P_G_lankle_cut  = P_G_lankle_cut (1:timeIndexToPeak2_f,:);
        P_G_lhip_cut    = P_G_lhip_cut   (1:timeIndexToPeak2_f,:);
        P_G_lsho_cut    = P_G_lsho_cut   (1:timeIndexToPeak2_f,:);
                   
        P_G_rtoe_cut    = P_G_rtoe_cut   (1:timeIndexToPeak2_f,:);
        P_G_rhee_cut    = P_G_rhee_cut   (1:timeIndexToPeak2_f,:);
        P_G_rankle_cut  = P_G_rankle_cut (1:timeIndexToPeak2_f,:);
        P_G_rhip_cut    = P_G_rhip_cut   (1:timeIndexToPeak2_f,:);
        P_G_rsho_cut    = P_G_rsho_cut   (1:timeIndexToPeak2_f,:);
                    
        P_G_tors_cut    = P_G_tors_cut  (1:timeIndexToPeak2_f,:);
        P_G_imuA_cut    = P_G_imuA_cut  (1:timeIndexToPeak2_f,:);
        P_G_imuB_cut    = P_G_imuB_cut  (1:timeIndexToPeak2_f,:);
        P_G_imuC_cut    = P_G_imuC_cut  (1:timeIndexToPeak2_f,:);
        
        
        halfIndex = round(length(accCut)*0.5);
        [~,timeIndexToPeak2_accl] = max(abs(accCut(halfIndex:end,chosenA_ID)));
        timeIndexToPeak2_accl = halfIndex + timeIndexToPeak2_accl;
        
        accCut = accCut(1:timeIndexToPeak2_accl,:);
        t_cut_imu = t_cut_imu(1:timeIndexToPeak2_accl);
        omegaCut = omegaCut(1:timeIndexToPeak2_accl,:);
        
        q1 = q1(1:timeIndexToPeak2_f,:);
        q2 = q2(1:timeIndexToPeak2_f,:);
        dq1 = dq1(1:timeIndexToPeak2_f,:);
        dq2 = dq2(1:timeIndexToPeak2_f,:);
        ddq1 = ddq1(1:timeIndexToPeak2_f,:);
        ddq2 = ddq2(1:timeIndexToPeak2_f,:);
    
        
        if(plotSecondCutData)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;
            
            subplot(411);
            plot(t_cut_vicon,fCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Force [N] ','FontSize',12);
            axis tight;
            grid on;
            title(sprintf('Subject : %d, Trial : %d, second cut',subjectID,trialID));
            
            subplot(412);
            plot(t_cut_vicon,momCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Moment [Nmm] ','FontSize',12);
            axis tight;
            grid on;
      
            subplot(413);
            plot(t_cut_imu,accCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('LinAcc [m/sec^2]','FontSize',12);
            axis tight;
            grid on;
            
            subplot(414);
            plot(t_cut_imu,omegaCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('AngVel [rad/s]','FontSize',12);
            axis tight;
            grid on;
        end
        
        if (plotJointQuantitiesSecondCut)
         fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;
        
            subplot(311);
            plot1 = plot(t_cut_vicon,q1.*(180/pi),'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(t_cut_vicon,q2.*(180/pi),'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$q_1$','$q_2$','Location','northeast');
            %title('Joint Quantities','FontSize',15);
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Angle [deg]','FontSize',15);
            axis tight;
            grid on;  
            title(sprintf('Subject %d, Trial %d, Joint Quantities secondCut',subjectID,trialID));
        
            subplot(312);
            plot1 = plot(t_cut_vicon,(180/pi)*dq1,'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(t_cut_vicon,(180/pi)*dq2,'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$\dot q_{1}$','$\dot q_{2}$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Velocity [deg/s]','FontSize',15);
            axis tight;
            grid on;
        
            subplot(313);
            plot1 = plot(t_cut_vicon,(180/pi)*ddq1,'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(t_cut_vicon,(180/pi)*ddq2,'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$\ddot q_{1}$','$\ddot q_{2}$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Acceleration [deg/s^2]','FontSize',15);
            axis tight;
            grid on;   
        end
        
        % 4.3 cutting settlingTimeCutIndex samples after first peak and before second peak
        
        fCut = fCut(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        t_cut_vicon = t_cut_vicon(settlingCutTimeIndex:(end-settlingCutTimeIndex));
        
        momCut = momCut(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        
        P_G_ltoe_cut    = P_G_ltoe_cut   (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_lhee_cut    = P_G_lhee_cut   (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_lankle_cut  = P_G_lankle_cut (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_lhip_cut    = P_G_lhip_cut   (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_lsho_cut    = P_G_lsho_cut   (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
                   
        P_G_rtoe_cut    = P_G_rtoe_cut   (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_rhee_cut    = P_G_rhee_cut   (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_rankle_cut  = P_G_rankle_cut (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_rhip_cut    = P_G_rhip_cut   (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_rsho_cut    = P_G_rsho_cut   (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
                    
        P_G_tors_cut    = P_G_tors_cut  (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_imuA_cut    = P_G_imuA_cut  (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_imuB_cut    = P_G_imuB_cut  (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_imuC_cut    = P_G_imuC_cut  (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        
        
        accCut = accCut(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        t_cut_imu = t_cut_imu(settlingCutTimeIndex:(end-settlingCutTimeIndex));
        
        omegaCut = omegaCut(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        
        q1 = q1(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        q2 = q2(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        dq1 = dq1(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        dq2 = dq2(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        ddq1 = ddq1(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        ddq2 = ddq2(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        
        if(length(t_cut_vicon) - length(t_cut_imu)) ~= 0;
            
            t_cutIndex = min(length(t_cut_vicon),length(t_cut_imu));
            
            t_cut_vicon = t_cut_vicon(1:t_cutIndex);
            fCut = fCut(1:t_cutIndex,:);
            momCut = momCut(1:t_cutIndex,:);
            
            t_cut_imu = t_cut_imu(1:t_cutIndex);
            accCut = accCut(1:t_cutIndex,:);
            omegaCut = omegaCut(1:t_cutIndex,:);
            
            P_G_ltoe_cut    = P_G_ltoe_cut   (1:t_cutIndex,:);
            P_G_lhee_cut    = P_G_lhee_cut   (1:t_cutIndex,:);
            P_G_lankle_cut  = P_G_lankle_cut (1:t_cutIndex,:);
            P_G_lhip_cut    = P_G_lhip_cut   (1:t_cutIndex,:);
            P_G_lsho_cut    = P_G_lsho_cut   (1:t_cutIndex,:);
                   
            P_G_rtoe_cut    = P_G_rtoe_cut   (1:t_cutIndex,:);
            P_G_rhee_cut    = P_G_rhee_cut   (1:t_cutIndex,:);
            P_G_rankle_cut  = P_G_rankle_cut  (1:t_cutIndex,:);
            P_G_rhip_cut    = P_G_rhip_cut   (1:t_cutIndex,:);
            P_G_rsho_cut    = P_G_rsho_cut   (1:t_cutIndex,:);
                    
            P_G_tors_cut    = P_G_tors_cut  (1:t_cutIndex,:);
            P_G_imuA_cut    = P_G_imuA_cut  (1:t_cutIndex,:);
            P_G_imuB_cut    = P_G_imuB_cut  (1:t_cutIndex,:);
            P_G_imuC_cut    = P_G_imuC_cut  (1:t_cutIndex,:);

            q1 = q1(1:t_cutIndex,:);
            q2 = q2(1:t_cutIndex,:);
            dq1 = dq1(1:t_cutIndex,:);
            dq2 = dq2(1:t_cutIndex,:);
            ddq1 = ddq1(1:t_cutIndex,:);
            ddq2 = ddq2(1:t_cutIndex,:);
            
        end
         
        

        
        if(plotSettlingTimeCutData)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;
            
            subplot(411);
            plot(t_cut_vicon,fCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Force [N] ','FontSize',12);
            axis tight;
            grid on;
            title(sprintf('Subject : %d, Trial : %d, cutted initial and final secs',subjectID,trialID));
             
            subplot(412);
            plot(t_cut_vicon,momCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Moment [Nmm] ','FontSize',12);
            axis tight;
            grid on; 
            
            subplot(413);
            plot(t_cut_imu,accCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('LinAcc [m/sec^2]','FontSize',12);
            axis tight;
            grid on;
            
            subplot(414);
            plot(t_cut_imu,omegaCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('AngVel [rad/s]','FontSize',12);
            axis tight;
            grid on;
        end
        
        if(plotJointQuantitiesSettlingTimeCut)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;
        
            subplot(311);
            plot1 = plot(t_cut_vicon,q1.*(180/pi),'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(t_cut_vicon,q2.*(180/pi),'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$q_1$','$q_2$','Location','northeast');
            %title('Joint Quantities','FontSize',15);
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Angle [deg]','FontSize',15);
            axis tight;
            grid on;  
            title(sprintf('Subject %d, Trial %d, Joint Quantities SettlingTimeCut',subjectID,trialID));
        
            subplot(312);
            plot1 = plot(t_cut_vicon,(180/pi)*dq1,'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(t_cut_vicon,(180/pi)*dq2,'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$\dot q_{1}$','$\dot q_{2}$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Velocity [deg/s]','FontSize',15);
            axis tight;
            grid on;
        
            subplot(313);
            plot1 = plot(t_cut_vicon,(180/pi)*ddq1,'lineWidth',1.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(t_cut_vicon,(180/pi)*ddq2,'lineWidth',1.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$\ddot q_{1}$','$\ddot q_{2}$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',15);
            xlabel('Time [s]','FontSize',15);
            ylabel('Acceleration [deg/s^2]','FontSize',15);
            axis tight;
            grid on;
        end
        
        %% 5. filtering 
               
        % ====force plate data (filtered in freq)
        f_filt = frequencyFilterSignal(fCut, cutOffFreqForcePlate, 1/samplingTime);
        mom_filt = frequencyFilterSignal(momCut, cutOffFreqForcePlate, 1/samplingTime);
        
        % ====imu data (filtered in freq)
        accl_filt = frequencyFilterSignal (accCut, cutOffFreqIMU, 1/samplingTime);
        omega_filt = frequencyFilterSignal (omegaCut, cutOffFreqIMU, 1/samplingTime);
        
        
        if(plotFilteredData)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            
            grid on;
            subplot(411);
            plot(t_cut_vicon,f_filt, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Force [N] ','FontSize',12);
            axis tight;
            grid on; 
            title(sprintf('Subject : %d, Trial : %d, filtered data',subjectID,trialID));
            
            subplot(412);
            plot(t_cut_vicon,mom_filt, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Moment [Nmm] ','FontSize',12);
            axis tight;
            grid on;
            
            subplot(413);
            plot(t_cut_vicon,accl_filt, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('LinAcc [m/sec^2]','FontSize',12);
            axis tight;
            grid on;
            subplot(414);
            plot(t_cut_vicon,omega_filt, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('AngVel [rad/s]','FontSize',12);
            axis tight;
            grid on;
        end 
        
        %% 6. data storing
        
        synchronisedData(subjectID,trialID).samplingTime = samplingTime;
        synchronisedData(subjectID,trialID).dataTime = t_cut_vicon;
        synchronisedData(subjectID,trialID).aLin_imu_imu = accl_filt;
        synchronisedData(subjectID,trialID).omega_imu_imu = omega_filt; 
        synchronisedData(subjectID,trialID).imu = [accl_filt accl_filt];
        synchronisedData(subjectID,trialID).wrench_fp_fp = [mom_filt,f_filt];
        
        synchronisedData(subjectID,trialID).P_G_ltoe =P_G_ltoe;
        synchronisedData(subjectID,trialID).P_G_lhee =P_G_lhee;
        synchronisedData(subjectID,trialID).P_G_lankle =P_G_lankle;
        synchronisedData(subjectID,trialID).P_G_lhip =P_G_lhip;
        synchronisedData(subjectID,trialID).P_G_lsho =P_G_lsho; 
        
        synchronisedData(subjectID,trialID).P_G_rtoe =P_G_rtoe;
        synchronisedData(subjectID,trialID).P_G_rhee = P_G_rhee;
        synchronisedData(subjectID,trialID).P_G_rankle =P_G_rankle;
        synchronisedData(subjectID,trialID).P_G_rhip = P_G_rhip;
        synchronisedData(subjectID,trialID).P_G_rsho = P_G_rsho;
        
        synchronisedData(subjectID,trialID).P_G_tors =P_G_tors;
        synchronisedData(subjectID,trialID).P_G_imuA = P_G_imuA;
        synchronisedData(subjectID,trialID).P_G_imuB =P_G_imuB;
        synchronisedData(subjectID,trialID).P_G_imuC =P_G_imuC;
        
        synchronisedData(subjectID,trialID).P_G_1 = P_G_1;
        synchronisedData(subjectID,trialID).P_G_2 = P_G_2;
        synchronisedData(subjectID,trialID).P_G_3 = P_G_3;
        
        synchronisedData(subjectID,trialID).q = [q1 q2];
        synchronisedData(subjectID,trialID).dq = [dq1 dq2];
        synchronisedData(subjectID,trialID).ddq = [ddq1 ddq2];  
        
    end   
    fprintf('\n');
 end

%% save result 
save('./experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat','synchronisedData');
fprintf('---------\n\n');
fprintf('Finished synchronising data\n');