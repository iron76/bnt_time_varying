% synchroniseCaptureData
% Script to synchronise the captured IMU and VICON/force plate data. The 
% It looks for the spikes in beginning and end of both the sources 
% (corresponding to the little tiptoe gesture performed by subjects in the 
% start and conclusion of the bowing motion experiment). It:
% 1. filters data;
% 2. interpolates data to obatain the same samplingTime for all of them;
% 3. syncrhonises data:
%    3.1 re-aligns data on the first peak cutting data before it;
%    3.2 cuts data after the second peak;
%    3.3 cuts for both peaks 1sec of acquisition 
% 4. converts eventually data in ISU;
% 5. stores data into a file.mat


% Assumptions:
% - C3D processing is used to generate the Vicon/Force plate data in .mat
% files;
% - previous processing to generate the IMU data in .mat file;


clc; clear; close all;

%% testOptions
plotRawData = false;
plotFirstCutData = false;
plotSecondCutData = false;
plotSettlingTimeCutData = false;
plotFinalConvertedData = false;

%% load data sources
load('./experiments/humanFixedBase/data/VICONsaveDataGen16.mat');
load('./experiments/humanFixedBase/data/imuExtractedDataGen16.mat');
addpath('./experiments/humanFixedBase/helperFunctions/');

subjectIDList =1:12;
trialIDList = 1:4;

% Setting parameters
samplingTime = 1e-2;  %uniform sampling rate 
settlingCutTime = 0.5;  %1sec
settlingCutTimeIndex = ceil(settlingCutTime/samplingTime); %ensuring integer variable

forcePlateSamplingTime = 1e-3;
imuSamplingTime = 1e-2;
markersSamplingTime = 1e-2;


%% iterate through each trial 
for subjectID = subjectIDList
    for trialID = trialIDList
        
        %% 1.filtering --> TO DO
        
        %% 2. data interpolation
        
        % wrench raw data
        f_raw = subjectData(subjectID,trialID).analogsFOR;
        mom_raw = subjectData(subjectID,trialID).analogsMOM;
        
        t_f_raw = 0:forcePlateSamplingTime:(forcePlateSamplingTime*(size(f_raw,1) -1));
        t_f = 0:samplingTime:t_f_raw(end);
        f = interp1(t_f_raw,f_raw,t_f);
        mom = interp1(t_f_raw,mom_raw,t_f);
        
        % imu raw data
        accl_raw = imuData(subjectID,trialID).accln;
        omega_raw = imuData(subjectID,trialID).gyro;
          
        if(~isempty(find(diff(imuData(subjectID,trialID).t)<=0)))
            ttemp = imuData(subjectID,trialID).t;
            t_imu_raw = linspace(0,ttemp(end),length(ttemp))./1000;
        else
            t_imu_raw = imuData(subjectID,trialID).t./1000;
        end
        
        t_imu = 0:samplingTime:t_imu_raw(end);
        accl = interp1(t_imu_raw,accl_raw,t_imu);
        omega = interp1(t_imu_raw,omega_raw,t_imu);

        if(plotRawData)
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
            title(sprintf('Subject : %d, Trial : %d',subjectID,trialID));
            
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
        
        
        % markers raw data 
        t_markers_raw = 0:markersSamplingTime:(markersSamplingTime*size(subjectData(subjectID,trialID).markers.ltoe,1));
        t_markers_raw = t_markers_raw(1:end-1);
        
        P_G_ltoe = interp1( t_markers_raw,subjectData(subjectID,trialID).markers.ltoe,t_f);
        P_G_lhee = interp1( t_markers_raw,subjectData(subjectID,trialID).markers.lhee,t_f);
        P_G_lankle = interp1( t_markers_raw,subjectData(subjectID,trialID).markers.lankle,t_f);
        P_G_lhip = interp1( t_markers_raw,subjectData(subjectID,trialID).markers.lhip,t_f);
        P_G_lsho =interp1( t_markers_raw,subjectData(subjectID,trialID).markers.lsho,t_f); 
        
        P_G_rtoe = interp1(t_markers_raw,subjectData(subjectID,trialID).markers.rtoe,t_f);
        P_G_rhee = interp1( t_markers_raw,subjectData(subjectID,trialID).markers.rhee,t_f);
        P_G_rankle =interp1(t_markers_raw,subjectData(subjectID,trialID).markers.rankle,t_f);
        P_G_rhip = interp1( t_markers_raw,subjectData(subjectID,trialID).markers.rhip,t_f);
        P_G_rsho = interp1(t_markers_raw,subjectData(subjectID,trialID).markers.rsho,t_f); 
        
        P_G_tors =interp1( t_markers_raw,subjectData(subjectID,trialID).markers.tors,t_f); 
        P_G_imuA = interp1( t_markers_raw,subjectData(subjectID,trialID).markers.imuA,t_f);
        P_G_imuB = interp1( t_markers_raw,subjectData(subjectID,trialID).markers.imuB,t_f);
        P_G_imuC =interp1( t_markers_raw,subjectData(subjectID,trialID).markers.imuC,t_f);
       
        %% 3. synchronising data
        
        % 3.1 synchronising peaks of force with imu 
        
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
            title(sprintf('Subject : %d, Trial : %d',subjectID,trialID));
            
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

        P_G_ltoe_cut    = P_G_ltoe   (timeIndexToPeak1_f:end,:);
        P_G_lhee_cut    = P_G_lhee   (timeIndexToPeak1_f:end,:);
        P_G_lankle_cut  = P_G_lankle (timeIndexToPeak1_f:end,:);
        P_G_lhip_cut    = P_G_lhip   (timeIndexToPeak1_f:end,:);
        P_G_lsho_cut    = P_G_lsho   (timeIndexToPeak1_f:end,:);
                   
        P_G_rtoe_cut    = P_G_rtoe   (timeIndexToPeak1_f:end,:);
        P_G_rhee_cut    = P_G_rhee   (timeIndexToPeak1_f:end,:);
        P_G_rankle_cut  = P_G_rankle  (timeIndexToPeak1_f:end,:);
        P_G_rhip_cut    = P_G_rhip   (timeIndexToPeak1_f:end,:);
        P_G_rsho_cut    = P_G_rsho   (timeIndexToPeak1_f:end,:);
                    
        P_G_tors_cut    = P_G_tors  (timeIndexToPeak1_f:end,:);
        P_G_imuA_cut    = P_G_imuA  (timeIndexToPeak1_f:end,:);
        P_G_imuB_cut    = P_G_imuB  (timeIndexToPeak1_f:end,:);
        P_G_imuC_cut    = P_G_imuC  (timeIndexToPeak1_f:end,:);
        
        
        
        % 3.2 cutting data after last peak 
        
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
        P_G_rankle_cut   = P_G_rankle_cut  (1:timeIndexToPeak2_f,:);
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
            
            subplot(412);
            plot(t_cut_vicon,momCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Moment [Nmm] ','FontSize',12);
            axis tight;
            grid on;
            title(sprintf('Subject : %d, Trial : %d',subjectID,trialID));
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
        
        % 3.3 cutting settlingTimeCutIndex samples after first peak and before second peak
        
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
        P_G_rankle_cut  = P_G_rankle_cut  (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_rhip_cut    = P_G_rhip_cut   (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_rsho_cut    = P_G_rsho_cut   (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
                    
        P_G_tors_cut    = P_G_tors_cut  (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_imuA_cut    = P_G_imuA_cut  (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_imuB_cut    = P_G_imuB_cut  (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        P_G_imuC_cut    = P_G_imuC_cut  (settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        
        
        accCut = accCut(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        t_cut_imu = t_cut_imu(settlingCutTimeIndex:(end-settlingCutTimeIndex));
        
        omegaCut = omegaCut(settlingCutTimeIndex:(end-settlingCutTimeIndex),:);
        
        if (length(t_cut_vicon) - length(t_cut_imu)) ~= 0;
            
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
            title(sprintf('Subject : %d, Trial : %d',subjectID,trialID));     
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
        
        %% 4. data conversion
        
        % VICON 
        % acquired data in mm --> converted in m
        convConst = 1e-3;
        P_G_ltoe =  convConst * P_G_ltoe_cut;
        P_G_lhee = convConst * P_G_lhee_cut;
        P_G_lankle = convConst * P_G_lankle_cut;
        P_G_lhip = convConst * P_G_lhip_cut;
        P_G_lsho = convConst * P_G_lsho_cut; 
        
        P_G_rtoe = convConst * P_G_rtoe_cut;
        P_G_rhee = convConst *  P_G_rhee_cut;
        P_G_rankle =convConst * P_G_rankle_cut;
        P_G_rhip =convConst * P_G_rhip_cut;
        P_G_rsho =convConst * P_G_rsho_cut;
        
        P_G_tors =convConst * P_G_tors_cut;
        P_G_imuA = convConst * P_G_imuA_cut;
        P_G_imuB =convConst * P_G_imuB_cut;
        P_G_imuC =convConst * P_G_imuC_cut;
        
        % FORCE PLATE
        % acquired force datain N --> no conversion;
        % acquired moment data in Nmm --> converted in Nm;
        momCut =  1e-3* momCut;
        
        % IMU
        % acquired acc data in m/s^2 --> no conversion;
        % acquired omega data in rad/s^2 --> no conversion;
        
        
        if(plotFinalConvertedData)
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
            title(sprintf('Subject : %d, Trial : %d',subjectID,trialID));
            
            subplot(412);
            plot(t_cut_vicon,momCut, 'lineWidth',2.0);
            xlabel('Time [s]','FontSize',12);
            ylabel('Moment [Nm] ','FontSize',12);
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
        
        %% 5. data storing
        
        synchronisedData(subjectID,trialID).samplingTime = samplingTime;
        synchronisedData(subjectID,trialID).dataTime = t_cut_vicon;
        synchronisedData(subjectID,trialID).aLin_imu_imu = accCut;
        synchronisedData(subjectID,trialID).omega_imu_imu = omegaCut; 
        synchronisedData(subjectID,trialID).wrench_fp_fp = [momCut,fCut];
        
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
        
    end
end

%% save result 
save('./experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat','synchronisedData');
