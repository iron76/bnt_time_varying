% load data

load('VICONsaveData.mat');
load('imuExtractedData.mat');

subjectIDList = [1];
trialIDList = 1:3;


for subjectID=1:length(subjectIDList)
    for trialID=1:length(trialIDList)
        
        %% compute VICON, IMU time difference
        figure;
        subplot(3,1,1);
        f_raw = subjectData(subjectID,trialID).grwsFOR;
        t_raw_vicon_f = 0:(1/100):(0.01*size(f_raw,1));
        t_raw_vicon_f = t_raw_vicon_f(1:end-1)./10;

        t_vicon = 0:(1/1000):t_raw_vicon_f(end);
        f = interp1(t_raw_vicon_f,f_raw,t_vicon);

        plot(t_vicon,f);axis tight;

        title(sprintf('Subject : %d, Trial : %d',subjectID,trialID));
        accl_raw = imuData(subjectID,trialID).accln;
        omega_raw = imuData(subjectID,trialID).gyro;
        
        if(~isempty(find(diff(imuData(subjectID,trialID).t)<=0)))
            ttemp = imuData(subjectID,trialID).t;
            t_imu_raw = linspace(0,ttemp(end),length(ttemp))./1000;
        else
            t_imu_raw = imuData(subjectID,trialID).t./1000;
        end
        
        

        t_imu = 0:(1/1000):t_imu_raw(end);
        accl = interp1(t_imu_raw,accl_raw,t_imu);
        omega = interp1(t_imu_raw,omega_raw,t_imu);

        subplot(3,1,2);plot(t_imu,accl);axis tight;

        [val,chosenF_ID] = max(f(1,:));
        [val,timeIndexToPeak1_f] = max(f(1:round(end*0.5),chosenF_ID));
        timeToPeak1_f = t_vicon(timeIndexToPeak1_f);

        [val,chosenA_ID] = max(accl(1,:));
        [val,timeIndexToPeak1_accl] = max(accl(1:round(end*0.5),chosenA_ID));

        timeToPeak1_accl = t_imu(timeIndexToPeak1_accl);

        indicesToShift = timeIndexToPeak1_accl - timeIndexToPeak1_f;
        delta_t_time = t_imu(indicesToShift);

        t_imu_shifted = t_imu(indicesToShift:end) - t_imu(indicesToShift);
        accl_shifted = accl(indicesToShift:end,:);
        omega_shifted = omega(indicesToShift:end,:);
        %interp1(t_imu,accl,t_imu_shifted);

        subplot(3,1,3);
        plot(t_imu_shifted,accl_shifted);
        axis tight;
        
        
        %% acceleration timeshift resultStore
        imu_vicon_shiftedData(subjectID,trialID).t_imu = t_imu_shifted;
        imu_vicon_shiftedData(subjectID,trialID).a_imu_imulin = accl_shifted;
        imu_vicon_shiftedData(subjectID,trialID).v_imu_imurot = omega_shifted;

        %% VICON interpolation, variable rename and resultStore
        
        t_raw_vicon_p = 0:(1/100):(0.01*size(subjectData(subjectID,trialID).markers.ltoe,1));
        t_raw_vicon_p = t_raw_vicon_p(1:end-1);
        
        imu_vicon_shiftedData(subjectID,trialID).t_vicon = t_vicon;
        
        imu_vicon_shiftedData(subjectID,trialID).P_G_ltoe = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.ltoe,t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_G_lhee = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.lhee,t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_G_lankle = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.lankle,t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_G_lhip = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.lhip,t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_G_lsho = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.lsho,t_vicon); 
        
        imu_vicon_shiftedData(subjectID,trialID).P_G_rtoe = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.rtoe,t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_G_rhee = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.rhee,t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_G_rankle = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.rankle,t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_G_rhip = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.rhip,t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_G_rsho = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.rsho,t_vicon); 
        
        imu_vicon_shiftedData(subjectID,trialID).P_G_tors = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.tors,t_vicon); 
        imu_vicon_shiftedData(subjectID,trialID).P_G_imuA = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.imuA,t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_G_imuB = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.imuB,t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_G_imuC = interp1( t_raw_vicon_p,subjectData(subjectID,trialID).markers.imuC,t_vicon);
        
        imu_vicon_shiftedData(subjectID,trialID).fx_PWAPWA_PWA = interp1( t_raw_vicon_f,[subjectData(subjectID,trialID).grwsMOM,subjectData(subjectID,trialID).grwsFOR],t_vicon);
        %imu_vicon_shiftedData(subjectID,trialID).f_GPWA_s = interp1( t_raw_vicon_f,[subjectData(subjectID,trialID).analogsFOR,subjectData(subjectID,trialID).analogsMOM],t_vicon);
        imu_vicon_shiftedData(subjectID,trialID).P_PWA_C = interp1( t_raw_vicon_f,subjectData(subjectID,trialID).grwsPOS,t_vicon);
        % extract points needed in correct variable names
      
    end
end

save('./experiments/humanFixedBase/IMU_VICON_ShiftedData.mat','imu_vicon_shiftedData');