
%ORGANISEBERDYCOMPATIBLESENSORDATA Reorganises the sensor data from human
%capture into a BERDY compatible form
%   loads sensor transfors and organises data to generate the y and ys
%   matrices in correct form

clc; clear; close all;

%% load processed data and sensor link transforms
load('./experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat');
load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');

%% selected subjects and trials
subjectList = 1;
trialList = 1;  

%% iterate through each computing transforms each time

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d ',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d ',trialID);

        currentTrial = synchronisedData(subjectID,trialID);
        
        data.dataTime = currentTrial.dataTime;
        data.q = currentTrial.q;
        data.dq = currentTrial.dq;
        data.ddq = currentTrial.ddq;
        
        
        %% data expressed in sensor frame --> data.y_sensFrame --> for Ymatrix created manually

        %IMU
        aLin_imu  = currentTrial.aLin_imu_imu;                %linear part of a, in imu frame
        %a_imu_imu   = [zeros(size(aLin_imu)), aLin_imu];   %twist in imu frame, ang-lin notation, assumption: aAng =0
        %data.ys_sensFrame_imu = a_imu_imu;

        %GYRO --> Assumption: considering that gyro can give us omegaDot
        omega_imu_imu  = currentTrial.omega_imu_imu;
        [omegaDot_imuX,~] = SgolayDerivation(3,57,omega_imu_imu(:,1),1e-2);
        [omegaDot_imuY,~] = SgolayDerivation(3,57,omega_imu_imu(:,2),1e-2);
        [omegaDot_imuZ,~] = SgolayDerivation(3,57,omega_imu_imu(:,3),1e-2);
        omegaDot_imu = [omegaDot_imuX omegaDot_imuY omegaDot_imuZ];
        
        a_imu_imu   = [omegaDot_imu, aLin_imu];
        data.ys_sensFrame_imu = a_imu_imu;
        
        %Force plate
        y_fp_fp = currentTrial.wrench_fp_fp;                %wrench in forceplate frame, ang-lin notation
        data.ys_sensFrame_fts = y_fp_fp;

        %% data expressed in link frame --> data.y_linkFrame
        
        currentTrialSens = sensorLinkTransforms(subjectID,trialID);
        
        % ====IMU in link2
        X_imu_2 = currentTrialSens.X_imu_2;
        X_2_imu = InverseAdjTransform(X_imu_2);
        a_2_imu = zeros(size(a_imu_imu))';
        
        for i = 1:length(data.dataTime)                    %twist in frame associate to link2
            a_2_imu(:,i)= X_2_imu * a_imu_imu(i,:)';
        end
                 
        data.ys_linkFrame_imu = a_2_imu';


        % ====Force plate in link0
        XStar_fp_0 = currentTrialSens.XStar_fp_0;
        XStar_0_fp = InverseAdjTransform(XStar_fp_0);
        y_0_fp = zeros(size(y_fp_fp))';
        
        for i = 1:length(data.dataTime)                    %wrench in frame associate to link0
            y_0_fp(:,i)= XStar_0_fp * y_fp_fp(i,:)';
        end
        
        data.ys_link0Frame_fts = y_0_fp';

        
        % ====Force plate in link1
        XStar_0_1 = currentTrialSens.XStar_0_1;
        y_1_fp = zeros (6,length(data.dataTime));
  
        for i = 1:length(data.dataTime)                    %wrench in frame associate to link1
            XStar_1_0{i} = InverseAdjTransform(XStar_0_1{i});
            y_1_fp(:,i)= XStar_1_0{i} * y_0_fp(:,i);
        end
        
        data.ys_link1Frame_fts = y_1_fp';
        
        
        %% Organising into a structure    
        BERDYFormattedSensorData(subjectID,trialID).data = data;       
    end
    fprintf('\n');
end

%% storing results
 save('./experiments/humanFixedBase/intermediateDataFiles/BERDYFormattedSensorData.mat','BERDYFormattedSensorData')
 
fprintf('---------\n');
fprintf('Done!\n');
