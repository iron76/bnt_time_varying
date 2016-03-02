%%COMPUTELINKSENSORFRAMES
% Script to compute the link-sensor frames associated with the
% Human-Dynamics estimation experiment. Sensors are an IMU
% placed on the chest and force place on the bottom of the foot. 
% The transforms being calculated are imu_X_2, and fp_XStar_2.
%
%
%% Note on the transforms : 
% The Featherstone convention is utilised [Featherstone(2008), 
% Rigid Body Dynamics Algorithms (2008), pg22]
% This implies that it is a transform for a spatial vector organised in 
% the angular-linear format.
%
%
%% Assumption of model URDF: 
% 1) P0, point in link0, origin of frame associated to the link0;
% 2) P1, point between link0 and link1, origin of frame associated to the link1.
% The orientation of P1 is the same of link1.
% 3) P2, point between link1 and link2, origin of frame associated to the link2.
% The orientation of P2 is the same of link2.
% 4) P3, point at the top of link2, since link3 doesn't exist.


clc; clear; close all;

%% load processed data and models
load('./experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat');
%isTest = 'false'; 

load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');

%% selected subjects and trials
subjectList = 1:12;
trialList = 1:4 ;  

%% iterate through each computing transforms each time

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d ',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d ',trialID);
        
        currentTrial = processedSensorData(subjectID,trialID);
        
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
        
        P_G_1 = currentTrial.P_G_1;
        P_G_2 = currentTrial.P_G_2;
        P_G_3 = currentTrial.P_G_3;

        [R_0_G,P_G_0] = computeFootRotation(P_G_lhee,P_G_rhee,P_G_ltoe,P_G_rtoe);    
        P_G_G = [0,0,0];
        P_fp_fp = [0,0,0];
         
        q1 = currentTrial.q1;
        q2 = currentTrial.q2;
        q = [q1 q2];
        
        %% Computing adjoint transform 2_X_imu  --> for Ymatrix
        % we need 2_X_imu for Y matrix . 
        % We want to compute imu_X_2 = imu_X_G * G_X_0 * 0_X_2
         
        currentModel = humanThreeLinkModelFromURDF(subjectID).dmodel;
        dmodel = currentModel;
        
        % Computing G_X_0 (const)
        r_G_from0toG = P_G_G - P_G_0 ;
        r_0_from0toG = R_0_G * r_G_from0toG';
        
        X_G_0 = computeAdjointTransform(R_0_G',r_0_from0toG);     
               
        % Computing imu_X_G (variant during motion). Since imu_X_G is used
        % to compute 2_X_imu that is const --> we are going to consider
        % only the mean of 10 samples for this computation.
        samples = 10;
       
        [R_G_imu, P_G_imu] = computeInitialIMURotation(P_G_imuA,P_G_imuB,P_G_imuC);
        r_G_fromGtoimu = P_G_imu - P_G_G;
        
        X_imu_G = computeAdjointTransform(R_G_imu',r_G_fromGtoimu); 
        
        % Computing 0_X_2     
        X_0_2 = AdjTransfFromLinkToRoot (dmodel, mean(q(1:samples,:)), 2);
        
        % Computing imu_X_2 
        X_imu_2 = X_imu_G * X_G_0 * X_0_2;

        
   
       %% notes: 
       % if we compute 0_X_1 and 1_X_2 as follows:
       %            
       % X_0_1 = AdjTransfFromLinkToRoot (humanThreeLinkModelFromURDF(subjectID).dmodel, mean(q(1:samples,:)), 1);                   
       % X_1_2   = inv(X_0_1) * X_0_2;
       %
       % we denote that 0_R_1 in 0_X_1 is not similar to an identity matrix
       % implying that a rotation has occurred (clearly the same condition is in 0_R_2).
       % Instead looking 1_R_2 in 1_X_2 is very similar to an identity matrix since there is not a
       % frame rotation.  
       % This is a consequence of Drake parsing.
       
       
        %% Computing adjoint transform fp_Xstar_0  --> for Ymatrix
        % we need fp_Xstar_0 for Y matrix
        
        %fixed rotation from Global and Force plate frames
        R_fp_G = [-1 0  0;
                   0 1  0;
                   0 0 -1];    
   
        %origin of fp frame in G frame (consider the force plate heigt -0.04330 m)       
        P_G_fp =  [0.23175,0.25400,-0.04330]; 
        
        %we want to compute r_fp_from0tofp
        r_G_from0tofp = P_G_fp - P_G_0;
        r_0_from0tofp = R_0_G * r_G_from0tofp';
        
        R_fp_0 = R_fp_G * R_0_G';
 
        XStar_fp_0 = [  R_fp_0    -R_fp_0*skew(r_0_from0tofp);
                       zeros(3)                R_fp_0         ];
        
           
        %% Computing adjoint transform 0_XStar_1  --> for Ymatrix
   
        XStar_0_1 = cell(size(q,1),1);
        
        for i = 1:size(q,1)
            XStar_0_1{i} = AdjTransfStarfFromLinkToRoot(dmodel, q(i,:),1);
        end
                  
        %% Organising into a structure    
        sensorLinkTransforms(subjectID,trialID).X_G_0 = X_G_0;
        sensorLinkTransforms(subjectID,trialID).X_imu_2 = X_imu_2;
        sensorLinkTransforms(subjectID,trialID).XStar_fp_0 = XStar_fp_0;
        sensorLinkTransforms(subjectID,trialID).XStar_0_1 = XStar_0_1; 
    end
    fprintf('\n');
end

%% storing results
%5if(strcmp(isTest,'true')~=1)
    save('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat','sensorLinkTransforms');
%end