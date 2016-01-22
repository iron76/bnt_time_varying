% computeLinkSensorFrames
% Script to compute the link-sensor frames associated with the
% Human-Dynamics estimation experiment. The 2 sensors are an IMU
% placed on the chest and force place on the bottom of my foot. 
% The transforms being calculated are 2_X_imu, and 0_X_fp and
% their corresponding inverses.
%
% Author: Naveen Kuppuswamy (naveen.kuppuswamy@iit.it)
% iCub Facility, Istituto Italiano di Tecnologia, 21 January 2016

%% load the synchronised dataset (including VICON and IMU data)
load('./experiments/humanFixedBase/intermediateDataFiles/synchronisedSensorData.mat');
isTest = 'false';

%% selected subjects and trials
subjectList = 1;
trialList = 1 ; 

%% iterate through each computing transforms each time
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
            
        P_G_1 = computeCentroidOfPoints(P_G_lankle,P_G_rankle);
        P_G_2 = computeCentroidOfPoints(P_G_lhip,P_G_rhip);
        P_G_3 = computeCentroidOfTriangle(P_G_lsho,P_G_rsho,P_G_tors);
        
        %% 0 frame
        [R_0_G,P_G_0] = computeFootRotation(P_G_lhee,P_G_rhee,P_G_ltoe,P_G_rtoe); 
       
        R_G_0 = R_0_G';
        R_G_1 = R_G_0;     % because point P_G_1 is fixed on the foot in URDF
        R_1_G = R_G_1';
        
        %% computing joint angles q1 and q2
                %% Computing q1 and q2  angles 
       
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

        %% Computing rotation matrix R_2_imuini
        
        R_0_1ini = euler2dcm([0,mean(q1(1:10)),0]); 
        R_1_2ini = euler2dcm([0,mean(q2(1:10)),0]); 
        R_G_2ini = R_G_0 * R_0_1ini * R_1_2ini;
        
        [R_G_imuini,P_G_Gimuini] = computeInitialIMURotation(P_G_imuA,P_G_imuB,P_G_imuC);
        R_2_imuini = R_G_2ini'* R_G_imuini;
        
        P_2_2imuini = P_G_Gimuini - mean(P_G_2(1:10,:));
        
        T_2_imu = [R_2_imuini, P_2_2imuini' .*1e-3 ; zeros(1,3) 1];
        X_2_imu = [R_2_imuini,                          zeros(3); ...
                  skew(P_2_2imuini' .*1e-3)*R_2_imuini, R_2_imuini];
        R_imuini_2 = R_2_imuini';      
        X_imu_2 = [R_imuini_2,                             zeros(3); ...
                   -R_2_imuini*skew(P_2_2imuini'.*1e-3),   R_imuini_2];   
      

        %% Xstar_0_fp
        %fixed rotation from Global and Force plate reference frames
        R_G_fp = [-1 0  0;
                   0 1  0; 
                   0 0 -1];   
              
        R_0_fp = R_0_G * R_G_fp;
        R_fp_0 = R_0_fp';
        
        % center of force plate in mm (below the force plate) in Global frame
        P_G_fp = [231.75,254,-43.3]; 
        
        r_G_from0toFp = P_G_0 - P_G_fp;
        r_0_from0toFp = R_0_G*r_G_from0toFp';
        
        r_0_from0toFpm = r_0_from0toFp*1e-3; %converting to m

       % XStar_0_fp = [R_0_fp' skew(r_0_from0toFpm)*R_0_fp'; zeros(3) R_0_fp'];
        XStar_0_fp = [R_0_fp skew(r_0_from0toFpm)*R_0_fp; zeros(3) R_0_fp];
        XStar_fp_0 = [R_fp_0,       -R_fp_0*skew(r_0_from0toFpm);
                      zeros(3),     R_fp_0];     
        
        R_1_0 = R_0_1ini';
        R_1_fp = R_1_0 * R_0_fp;
        r_G_from1toFpm = mean(P_G_1(1:10,:)) - P_G_fp;
        R_1_G = R_G_1';
        r_1_from1toFpm = R_1_G* r_G_from1toFpm';
        XStar_1_fp = [R_1_fp    skew(r_1_from1toFpm*1e-3)*R_1_fp;...
                     zeros(3)   R_1_fp];
        R_fp_1 = R_1_fp';
        XStar_fp_1 = [R_fp_1    -R_fp_1*skew(r_1_from1toFpm*1e-3);...
                      zeros(3)  R_fp_1];

        %% Organising into a structure          
        sensorLinkTransforms(subjectID,trialID).X_2_imu = X_2_imu;
        sensorLinkTransforms(subjectID,trialID).X_imu_2 = X_imu_2;
        sensorLinkTransforms(subjectID,trialID).XStar_0_fp = XStar_0_fp;
        sensorLinkTransforms(subjectID,trialID).XStar_fp_0 = XStar_fp_0;
        sensorLinkTransforms(subjectID,trialID).XStar_1_fp = XStar_1_fp;
        sensorLinkTransforms(subjectID,trialID).XStar_fp_1 = XStar_fp_1;
    end
    fprintf('\n');
end

%% save result in the preProcessingDataFiles folder
if(strcmp(isTest,'true')~=1)
    save('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat','sensorLinkTransforms');
end