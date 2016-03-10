
% TO DO: SCRIPT DESCRIPTION
%this analysis must to be computed only on one trial!

clc; clear; close all;

%% load data
load('./experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat');

subjectList = 1:12;
trialList = 1 ; %fixed trial!

for subjectID = subjectList
    for trialID = trialList
       
        currentTrial = processedSensorData(subjectID,trialID);
        sSelec = 30; %selected first 30 samples for the mean
        
        P_G_lheeM = mean(currentTrial.P_G_lhee(1:sSelec,:));
        P_G_ltoeM = mean(currentTrial.P_G_ltoe(1:sSelec,:));
        P_G_rheeM = mean(currentTrial.P_G_rhee(1:sSelec,:));
        P_G_rtoeM = mean(currentTrial.P_G_rtoe(1:sSelec,:));
        P_G_lankleM = mean(currentTrial.P_G_lankle(1:sSelec,:));
        P_G_lhipM = mean(currentTrial.P_G_lhip(1:sSelec,:));
        P_G_lshoM = mean(currentTrial.P_G_lsho(1:sSelec,:));
        P_G_rankleM = mean(currentTrial.P_G_rankle(1:sSelec,:));
        P_G_rhipM = mean(currentTrial.P_G_rhip(1:sSelec,:));
        P_G_rshoM = mean(currentTrial.P_G_rsho(1:sSelec,:));
        P_G_torsM = mean(currentTrial.P_G_tors(1:sSelec,:));
        P_G_imuAM = mean(currentTrial.P_G_imuA(1:sSelec,:));
        P_G_imuBM = mean(currentTrial.P_G_imuB(1:sSelec,:));
        P_G_imuCM = mean(currentTrial.P_G_imuC(1:sSelec,:));
        
        %% computing points P_G_1,P_G_2,P_G_3
        
        P_G_1M = computeCentroidOfPoints(P_G_lankleM,P_G_rankleM);
        P_G_2M = computeCentroidOfPoints(P_G_lhipM,P_G_rhipM);
        P_G_3M= computeCentroidOfPoints(P_G_lshoM,P_G_rshoM);
        
        %% computing alfa,beta,gamma sizes for URDF bounding boxes                
       
        % FOOT BOX SIZE
        footAlpha = max([norm(P_G_lheeM - P_G_rheeM),norm(P_G_ltoeM - P_G_rtoeM)]);      
        footBeta = P_G_1M(:,3); %height of P_G_1
        footGamma = max([norm(P_G_lheeM - P_G_ltoeM),norm(P_G_rheeM - P_G_rtoeM)]); 
        
        % LEG BOX SIZE
        legAlpha = norm(P_G_lankleM - P_G_rankleM);     
        legBeta =  P_G_2M(:,3) - P_G_1M(:,3);
        legGamma = 0.5*footGamma;
       
        % TORSO BOX SIZE
        torsoAlpha = norm(P_G_lshoM - P_G_rshoM);
        torsoBeta = P_G_3M(:,3) - P_G_2M(:,3);
        torsoGamma = legGamma;
         
        %% computing box volume and mass
        
        footV = footAlpha*footBeta*footGamma; 
        legV = legAlpha*legBeta*legGamma;
        torsoV = torsoAlpha*torsoBeta*torsoGamma;
        totalV = footV+legV+torsoV;
        
        f_fp_fp = currentTrial.wrench_fp_fp(:,4:6);
        f_fp_fpM = f_fp_fp(1:sSelec,:);
        
        totalMass = abs(mean(f_fp_fpM(:,3)) ./ 9.81);
        footMass = totalMass * (footV/totalV);
        legMass= totalMass * (legV/totalV);
        torsoMass = totalMass * (torsoV/totalV);
        
        %% computing inertias according to the URDF reference frames
       
        % inertias are computed considering homogeneous density of boxes
        footIxx = (footMass/12)*(footAlpha^2 + footBeta^2);
        footIyy = (footMass/12)*(footBeta^2 + footGamma^2);
        footIzz = (footMass/12)*(footGamma^2 + footAlpha^2);
        
        legIxx = (legMass/12)*(legAlpha^2 + legBeta^2);
        legIyy = (legMass/12)*(legBeta^2 + legGamma^2);
        legIzz = (legMass/12)*(legGamma^2 + legAlpha^2);

        torsoIxx = (torsoMass/12)*(torsoAlpha^2 + torsoBeta^2);
        torsoIyy = (torsoMass/12)*(torsoBeta^2 + torsoGamma^2);
        torsoIzz = (torsoMass/12)*(torsoGamma^2 + torsoAlpha^2);


        fprintf('Subject: %d \n',subjectID);
        fprintf('FOOT   Alpha: %3.3f,    Beta: %3.3f,    Gamma: %3.3f\n',footAlpha,footBeta,footGamma)
        fprintf('LEG    Alpha: %3.3f,    Beta: %3.3f,    Gamma: %3.3f\n',legAlpha,legBeta,legGamma);
        fprintf('TORSO  Alpha: %3.3f,    Beta: %3.3f,    Gamma: %3.3f\n',torsoAlpha,torsoBeta,torsoGamma);
        fprintf('MASS   Total: %3.3f,   Foot : %3.3f,    Leg : %3.3f,    Torso : %3.3f\n',totalMass,footMass,legMass,torsoMass);
        fprintf('FootINERTIAS   Ixx: %3.3f,  Iyy: %3.3f,  Izz: %3.3f\n',footIxx,footIyy,footIzz);
        fprintf('LegINERTIAS    Ixx: %3.3f,  Iyy: %3.3f,  Izz: %3.3f\n',legIxx,legIyy,legIzz);
        fprintf('TorsoINERTIAS  Ixx: %3.3f,  Iyy: %3.3f,  Izz: %3.3f\n',torsoIxx,torsoIyy,torsoIzz);
        fprintf('---------------------------------\n');
          
        %% data storing
        
        subjectParams(subjectID,trialID).footWidth = footAlpha;
        subjectParams(subjectID,trialID).footHeight = footBeta;
        subjectParams(subjectID,trialID).footDepth = footGamma;
        subjectParams(subjectID,trialID).legWidth = legAlpha;
        subjectParams(subjectID,trialID).legHeight = legBeta;
        subjectParams(subjectID,trialID).legDepth = legGamma;
        subjectParams(subjectID,trialID).torsoWidth = torsoAlpha;
        subjectParams(subjectID,trialID).torsoHeight = torsoBeta;
        subjectParams(subjectID,trialID).torsoDepth = torsoGamma;
      
        subjectParams(subjectID,trialID).footMass = footMass;
        subjectParams(subjectID,trialID).footIxx = footIxx;
        subjectParams(subjectID,trialID).footIyy = footIyy; 
        subjectParams(subjectID,trialID).footIzz = footIzz;
        
        subjectParams(subjectID,trialID).torsoMass = torsoMass;
        subjectParams(subjectID,trialID).torsoIxx = torsoIxx;
        subjectParams(subjectID,trialID).torsoIyy = torsoIyy;
        subjectParams(subjectID,trialID).torsoIzz = torsoIzz;
        
        subjectParams(subjectID,trialID).legMass = legMass;
        subjectParams(subjectID,trialID).legIxx = legIxx;
        subjectParams(subjectID,trialID).legIyy = legIyy;
        subjectParams(subjectID,trialID).legIzz = legIzz;

    end
end

save('./experiments/humanFixedBase/data/subjectSizeParams.mat','subjectParams');




