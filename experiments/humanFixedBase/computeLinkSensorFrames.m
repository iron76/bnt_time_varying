% computeLinkSensorFrames
% Script to compute the link-sensor frames associated with the
% Human-Dynamics estimation experiment. The 2 sensors are an IMU
% placed on the chest and force place on the bottom of my foot. 
% The transforms being calculated are 2_X_imu, and 0_X_fp and
% their corresponding inverses.
%
% %% Note on the transforms : 
% The Featherstone convention is utilised [Featherstone(2008), 
% Rigid Body Dynamics Algorithms (2008), pg22]
% This implies that its a transform for a spatial vector organised in 
% the angular-linear format.
%
% The motion vector transform below follows the formula given below
% X_B_A = [R_B_A  0; (-R_B_A)(px) R_B_A]
%
% where if A is located at O and B at P and p is the
% position vector OP in A coordinates, and R_B_A is the 
% rotation matrix transforming from A to B coordinates    
% Corresponding invere is given by,
% X_A_B = [R_B_A' 0; (px)(R_B_A') R_B_A']
%
% The force vector transform below follows the formula given below
% XStar_B_A = [R_B_A  (-R_B_A)(px); 0  R_B_A]  
%
% where if A is located at O and B at P and p is the
% position vector OP in A coordinates, and R_B_A is the 
% rotation matrix transforming from A to B coordinates    
% Corresponding inverse is given by,
% XStar_A_B = [R_B_A' (px)(R_B_A'); 0 R_B_A']
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
        [R_G_0,P_G_0] = computeFootRotation(P_G_lhee,P_G_rhee,P_G_ltoe,P_G_rtoe); 
       
        R_0_G = R_G_0';
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
        
        P_G_2imuini = P_G_Gimuini - mean(P_G_2(1:10,:));
        P_2_2imuini = R_G_2ini'*P_G_2imuini';
        
        
        
        R_imuini_2 = R_2_imuini';      

        X_imu_2 = [R_imuini_2,                             zeros(3); ...
                   -R_imuini_2*skew(P_2_2imuini.*1e-3),   R_imuini_2];   
        X_2_imu = [R_2_imuini,                          zeros(3); ...
                  skew(P_2_2imuini.*1e-3)*R_2_imuini, R_2_imuini];
        
       
        %% Xstar_0_fp
        %fixed rotation from Global and Force plate reference frames
        R_G_fp = [-1 0  0; 0 1  0; 0 0 -1]; 
               
        R_0_fp = R_0_G * R_G_fp;
        R_fp_0 = R_0_fp';
        
        % center of force plate in mm (below the force plate) in Global frame
        P_G_fp = [231.75,254,-43.3]; 
        
        R_1_0 = R_0_1ini';
        R_1_fp = R_1_0 * R_0_fp;
        r_G_from1toFpm = P_G_fp - mean(P_G_1(1:10,:)) ;
        R_1_G = R_G_1';
        r_1_from1toFpm = R_1_G* r_G_from1toFpm';
        XStar_1_fpini = [R_1_fp    skew(r_1_from1toFpm*1e-3)*R_1_fp;...
                     zeros(3)   R_1_fp];
        R_fp_1 = R_1_fp';
        XStar_fp_1ini = [R_fp_1    -R_fp_1*skew(r_1_from1toFpm*1e-3);...
                      zeros(3)  R_fp_1];

        %% time varying version of  Xstar_0_fp
        R_0_1t = cell(size(q1));
        R_1_0t = cell(size(q1));
        R_1_fpt = cell(size(q1));
        XStar_1_fpt = cell(size(q1));
        XStar_fp_1t = cell(size(q1));
        for i = 1: length(q1)
            R_0_1t{i} = euler2dcm([0,q1(i),0]); 
            R_1_0t{i} = R_0_1t{i}';
            R_1_fpt{i}= R_1_0t{i} * R_0_fp;
            
            XStar_1_fpt{i} = [R_1_fpt{i}    skew(r_1_from1toFpm*1e-3)*R_1_fpt{i};...
                            zeros(3)   R_1_fpt{i}];
            R_fp_1t = R_1_fpt{i}';
            XStar_fp_1t{i} = [R_fp_1t    -R_fp_1t*skew(r_1_from1toFpm*1e-3);...
                             zeros(3)  R_fp_1t];
        end
                  
                  %% Organising into a structure          
        sensorLinkTransforms(subjectID,trialID).X_2_imu = X_2_imu;
        sensorLinkTransforms(subjectID,trialID).X_imu_2 = X_imu_2;
        sensorLinkTransforms(subjectID,trialID).XStar_1_fpini = XStar_1_fpini;
        sensorLinkTransforms(subjectID,trialID).XStar_fp_1ini = XStar_fp_1ini;
        sensorLinkTransforms(subjectID,trialID).XStar_1_fpt = XStar_1_fpt;
        sensorLinkTransforms(subjectID,trialID).XStar_fp_1t = XStar_fp_1t;
    end
    fprintf('\n');
end

%% save result in the preProcessingDataFiles folder
if(strcmp(isTest,'true')~=1)
    save('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat','sensorLinkTransforms');
end