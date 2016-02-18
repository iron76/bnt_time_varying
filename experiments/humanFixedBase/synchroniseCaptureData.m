% synchroniseCaptureData
% Script to synchronise the captured IMU and VICON data. 
% It is assumed that some kind of C3D processing is used to generate the 
% data in mat files. The script basically looks for the spikes in
% beginning and end of both the sources (corresponding to the little
% tiptoe gesture performed by subjects in the start and conclusion of
% the bowing motion experiment). The script then realigns the time indices
% to synchronise both the data sources and outputs a single matlab 
% structure that corresponds to the data.
%
% Author: Naveen Kuppuswamy (naveen.kuppuswamy@iit.it)
% iCub Facility, Istituto Italiano di Tecnologia, 21 January 2016

clear; 
close all;


%% load data sources
load('./experiments/humanFixedBase/data/newData/VICONsaveDataGen16.mat');
load('./experiments/humanFixedBase/data/newData/imuExtractedData.mat');
addpath('./experiments/humanFixedBase/helperFunctions/');

subjectIDList = 1;
trialIDList = 1;%1:4;

%% Sampling and filter settings
samplingFrequency = (1/100);

filterOpt.imu.toFilter = 'yes';
filterOpt.imu.window = 57;
filterOpt.imu.order = 3;

filterOpt.vicon.toFilter = 'no';
filterOpt.vicon.window = 23;
filterOpt.vicon.order = 4;

filterOpt.forcePlate.toFilter = 'yes';
filterOpt.forcePlate.window = 57;
filterOpt.forcePlate.order = 3;


%% iterate through each trial c
for subjectID = subjectIDList
    for trialID = trialIDList
        
        %% Prefiltering
        [viconDataFilt,imuDataFilt] = filterData(subjectData(subjectID,trialID),imuData(subjectID,trialID),filterOpt);

        %% compute VICON, IMU time difference
        figure;
        subplot(2,2,1);
        
        f_raw = viconDataFilt.analogsFOR;
        
        t_raw_vicon_f = 0:(1/1000):((1/1000)*size(f_raw,1));
        t_raw_vicon_f = t_raw_vicon_f(1:end-1);
        t_vicon = 0:(samplingFrequency):t_raw_vicon_f(end);
        f = interp1(t_raw_vicon_f,f_raw,t_vicon);

        plot(t_vicon,f);
        xlabel('t (sec)','FontSize',15);
        ylabel('A','FontSize',15);
        axis tight;
        grid on;
      
        title(sprintf('Subject : %d, Trial : %d',subjectID,trialID));
        accl_raw = imuDataFilt.accln;
        omega_raw = imuDataFilt.gyro;
        
        if(~isempty(find(diff(imuDataFilt.t)<=0)))
            ttemp = imuData(subjectID,trialID).t;
            t_imu_raw = linspace(0,ttemp(end),length(ttemp))./1000;
        else
            t_imu_raw = imuDataFilt.t./1000;
        end
      
        t_imu = 0:(samplingFrequency):t_imu_raw(end);
        accl = interp1(t_imu_raw,accl_raw,t_imu);
        omega = interp1(t_imu_raw,omega_raw,t_imu);

        subplot(2,2,3);
        plot(t_imu,accl);
        axis tight;
        grid on;
        xlabel('t (sec)','FontSize',15);
        ylabel('B','FontSize',15);
       % legend('1','2','3');
        
        %[~,chosenF_ID] = max(f(1,:));
        chosenF_ID = 3;
        [~,timeIndexToPeak1_f] = max(f(1:round(end*0.5),chosenF_ID));
        timeToPeak1_f = t_vicon(timeIndexToPeak1_f);

        [~,chosenA_ID] = max(accl(1,:));
        [val,timeIndexToPeak1_accl] = max(accl(1:round(end*0.5),chosenA_ID));

        timeToPeak1_accl = t_imu(timeIndexToPeak1_accl);

        if(timeToPeak1_accl>timeToPeak1_f)
            indicesToShift = timeIndexToPeak1_accl - timeIndexToPeak1_f;
            delta_t_time = t_imu(indicesToShift);

            t_imu_shifted = t_imu(indicesToShift:end) - t_imu(indicesToShift);
            accl_shifted = accl(indicesToShift:end,:);
            omega_shifted = omega(indicesToShift:end,:);

            %% acceleration timeshift resultStore
 
            t_viconShifted = t_vicon;

        else
            indicesToShift = timeIndexToPeak1_accl - timeIndexToPeak1_f;
            indicesToShift = abs(indicesToShift);
            delta_f_time = t_vicon(indicesToShift);
            
            t_viconShifted = t_vicon(indicesToShift:end) - t_vicon(indicesToShift);
            t_vicon = t_vicon(indicesToShift:end);

            t_imu_shifted = t_imu;
            accl_shifted = accl;
            omega_shifted = omega;
            
            tSelLen_imu = length(t_imu_shifted);
            tSelLen_vicon = length(t_viconShifted);
            
            % incase lengths are not same (i.e. one has longer timeseries)
            % then truncate the other
            if(tSelLen_vicon>tSelLen_imu)
                t_vicon = t_vicon(1:tSelLen_imu);
                t_viconShifted = t_viconShifted(1:tSelLen_imu);
            else
                t_imu_shifted = t_imu_shifted(1:tSelLen_vicon);
                accl_shifted = accl(1:tSelLen_vicon,:);
                omega_shifted = omega(1:tSelLen_vicon,:);            
            end
            
        end
        subplot(2,2,4);
        plot(t_imu_shifted,accl_shifted);
        axis tight;
        grid on;
        ylabel('Ashifted','FontSize',15);
        
        

        %% VICON interpolation, variable rename and resultStore
        t_raw_vicon_p = 0:(1/100):(0.01*size(viconDataFilt(subjectID,trialID).markers.ltoe,1));
        t_raw_vicon_p = t_raw_vicon_p(1:end-1);

        
        synchronisedData(subjectID,trialID).t_vicon = t_viconShifted;
        synchronisedData(subjectID,trialID).t_imu = t_imu_shifted;
        synchronisedData(subjectID,trialID).a_imu_imulin = accl_shifted;
        synchronisedData(subjectID,trialID).v_imu_imurot = omega_shifted;
        
        %original data in mm --> converted in m
        synchronisedData(subjectID,trialID).P_G_ltoe =1e-3* interp1( t_raw_vicon_p,viconDataFilt.markers.ltoe,t_vicon);
        synchronisedData(subjectID,trialID).P_G_lhee = 1e-3*interp1( t_raw_vicon_p,viconDataFilt.markers.lhee,t_vicon);
        synchronisedData(subjectID,trialID).P_G_lankle =1e-3* interp1( t_raw_vicon_p,viconDataFilt.markers.lankle,t_vicon);
        synchronisedData(subjectID,trialID).P_G_lhip = 1e-3*interp1( t_raw_vicon_p,viconDataFilt.markers.lhip,t_vicon);
        synchronisedData(subjectID,trialID).P_G_lsho =1e-3* interp1( t_raw_vicon_p,viconDataFilt.markers.lsho,t_vicon); 
        
        synchronisedData(subjectID,trialID).P_G_rtoe =1e-3* interp1( t_raw_vicon_p,viconDataFilt.markers.rtoe,t_vicon);
        synchronisedData(subjectID,trialID).P_G_rhee = 1e-3*interp1( t_raw_vicon_p,viconDataFilt.markers.rhee,t_vicon);
        synchronisedData(subjectID,trialID).P_G_rankle =1e-3* interp1( t_raw_vicon_p,viconDataFilt.markers.rankle,t_vicon);
        synchronisedData(subjectID,trialID).P_G_rhip = 1e-3*interp1( t_raw_vicon_p,viconDataFilt.markers.rhip,t_vicon);
        synchronisedData(subjectID,trialID).P_G_rsho = 1e-3*interp1( t_raw_vicon_p,viconDataFilt.markers.rsho,t_vicon); 
        
        synchronisedData(subjectID,trialID).P_G_tors =1e-3* interp1( t_raw_vicon_p,viconDataFilt.markers.tors,t_vicon); 
        synchronisedData(subjectID,trialID).P_G_imuA = 1e-3*interp1( t_raw_vicon_p,viconDataFilt.markers.imuA,t_vicon);
        synchronisedData(subjectID,trialID).P_G_imuB = 1e-3*interp1( t_raw_vicon_p,viconDataFilt.markers.imuB,t_vicon);
        synchronisedData(subjectID,trialID).P_G_imuC =1e-3* interp1( t_raw_vicon_p,viconDataFilt.markers.imuC,t_vicon);
        
        %% multiplying moments alone with 1e-3 to take Nmm to Nm
        synchronisedData(subjectID,trialID).f_fp = interp1( t_raw_vicon_f,[viconDataFilt.analogsMOM .* 1e-3,viconDataFilt.analogsFOR],t_vicon);

        subplot(2,2,2);
        plot(t_vicon,synchronisedData(subjectID,trialID).f_fp(:,4:6));
        axis tight;
        grid on;
        ylabel('Bshifted','FontSize',15);
        
      
    end
end

%% save result in the preProcessingDataFiles folder
save('./experiments/humanFixedBase/intermediateDataFiles/synchronisedSensorData.mat','synchronisedData');