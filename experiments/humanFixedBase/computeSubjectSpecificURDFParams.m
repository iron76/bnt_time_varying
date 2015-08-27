% load data

load('IMU_VICON_ShiftedData.mat');


subjectList = 1:3;
trialList = 1 ; 

for subjectID = 1:length(subjectList)
    for trialID = 1:length(trialList)
       temp = imu_vicon_shiftedData(subjectID,trialID);
        
       pSelec = 25;
        
       P_G_ltoe = mean(temp.P_G_ltoe(1:pSelec,:) );
       P_G_lhee = mean(temp.P_G_lhee(1:pSelec,:) );
       P_G_lankle = mean(temp.P_G_lankle(1:pSelec,:));
       P_G_lhip = mean(temp.P_G_lhip(1:pSelec,:));
       P_G_lsho = mean(temp.P_G_lsho(1:pSelec,:));
        
       P_G_rtoe = mean(temp.P_G_rtoe(1:pSelec,:) );
       P_G_rhee = mean(temp.P_G_rhee(1:pSelec,:) );
       P_G_rankle = mean(temp.P_G_rankle(1:pSelec,:));
       P_G_rhip = mean(temp.P_G_rhip(1:pSelec,:));
       P_G_rsho = mean(temp.P_G_rsho(1:pSelec,:));
        
       P_G_tors = mean(temp.P_G_tors(1:pSelec,:));
        
       fx_PWAPWA_1 = mean(temp.fx_PWAPWA_1(1:pSelec,:));
                
       P_G_1 = computeCentroidOfPoints(P_G_lankle,P_G_rankle);
       P_G_2 = computeCentroidOfPoints(P_G_lhip,P_G_rhip);
       P_G_3 = computeCentroidOfTriangle(P_G_lsho,P_G_rsho,P_G_tors);
        
       
       a = max([norm(P_G_lankle - P_G_lhee),norm(P_G_rankle - P_G_rhee)]).*1e-3;
       b = max([norm(P_G_lhee - P_G_ltoe),norm(P_G_rhee - P_G_rtoe)]).*1e-3;
       c = max([norm(P_G_rankle - P_G_rtoe),norm(P_G_lankle - P_G_ltoe)]).*1e-3;
       Her = computeHeronFormula(a,b,c); 
       hGamma = 2*((Her)./b);
       %CLA: adding Heron formula because we want the height of the ankle
       %in respect to the foot plane
       
       % FOOT BOX SIZE 
       footAlpha = max([norm(P_G_lhee - P_G_rhee),norm(P_G_ltoe - P_G_rtoe)]).*1e-3;
               %footBeta =  max([norm(P_G_lankle - P_G_lhee),norm(P_G_rankle - P_G_rhee)]).*1e-3;
       footBeta = hGamma;
       footGamma = max([norm(P_G_lhee - P_G_ltoe),norm(P_G_rhee - P_G_rtoe)]).*1e-3;
       
        footV = footAlpha*footBeta*footGamma;
       
       % LEG BOX SIZE
               %legAlpha = 0.5*(norm(P_G_rankle - P_G_lankle) + norm(P_G_rhip - P_G_lhip)).*1e-3; 
       legAlpha = norm(P_G_lankle - P_G_rankle).*1e-3;        
       legBeta =  max([norm(P_G_lhip - P_G_lankle),norm(P_G_rhip - P_G_rankle)]).*1e-3;
       legGamma = 0.5*(0.5*(norm(P_G_rhee - P_G_rtoe) + norm(P_G_lhee - P_G_ltoe))).*1e-3;
       
        legV = legAlpha*legBeta*legGamma;
       
       % TORSO BOX SIZE
                %torsoAlpha = 0.5*(norm(P_G_lhip - P_G_rhip)+norm(P_G_lsho - P_G_rsho)).*1e-3;
       torsoAlpha = norm(P_G_lsho - P_G_rsho).*1e-3;
       torsoBeta =  max(norm((P_G_lsho - P_G_lhip)),norm((P_G_rsho - P_G_rhip))).*1e-3;
       torsoGamma = legGamma;
       
        torsoV = torsoAlpha*torsoBeta*torsoGamma;
       
        
        
       totalV = footV+legV+torsoV;
       Totalm = norm(fx_PWAPWA_1(:,1:3)) ./ 9.81;
       footm = Totalm * (footV/totalV);
       legm = Totalm * (legV/totalV);
       torsom = Totalm * (torsoV/totalV);
       
       
       % inertias are computed according to the URDF reference frames
        footIxx = (footm/12)*(footAlpha^2 + footBeta^2);
        footIyy = (footm/12)*(footBeta^2 + footGamma^2);
        footIzz = (footm/12)*(footGamma^2 + footAlpha^2);
        

        legIxx = (legm/12)*(legAlpha^2 + legBeta^2);
        legIyy = (legm/12)*(legBeta^2 + legGamma^2);
        legIzz = (legm/12)*(legGamma^2 + legAlpha^2);


        torsoIxx = (torsom/12)*(torsoAlpha^2 + torsoBeta^2);
        torsoIyy = (torsom/12)*(torsoBeta^2 + torsoGamma^2);
        torsoIzz = (torsom/12)*(torsoGamma^2 + torsoAlpha^2);


       fprintf('Subject ID : %d \n',subjectID);
       fprintf('Foot Alpha : %3.3f, Beta : %3.3f, Gamma : %3.3f,\n',footAlpha,footBeta,footGamma)
       fprintf('Leg Alpha : %3.3f, Beta : %3.3f, Gamma : %3.3f,\n',legAlpha,legBeta,legGamma);
       fprintf('Torso Alpha : %3.3f, Beta : %3.3f, Gamma : %3.3f,\n',torsoAlpha,torsoBeta,torsoGamma);
       fprintf('Mass Total : %3.3f, Foot : %3.3f, Leg : %3.3f, Torso : %3.3f,\n',Totalm,footm,legm,torsom);
       fprintf('FootInertias : Ixx : %3.3f, Iyy : %3.3f, Izz : %3.3f\n',footIxx,footIyy,footIzz);
       fprintf('LegInertias : Ixx : %3.3f, Iyy : %3.3f, Izz : %3.3f\n',legIxx,legIyy,legIzz);
       fprintf('TorsoInertias : Ixx : %3.3f, Iyy : %3.3f, Izz : %3.3f\n',torsoIxx,torsoIyy,torsoIzz);
       fprintf('---------------------------------\n');
       
       subjectParams(subjectID,trialID).footAlpha = footAlpha;
       subjectParams(subjectID,trialID).footBeta = footBeta;
       subjectParams(subjectID,trialID).footGamma = footGamma;
       subjectParams(subjectID,trialID).legAlpha = legAlpha;
       subjectParams(subjectID,trialID).legBeta = legBeta;
       subjectParams(subjectID,trialID).legGamma = legGamma;
       subjectParams(subjectID,trialID).torsoAlpha = torsoAlpha;
       subjectParams(subjectID,trialID).torsoBeta = torsoBeta;
       subjectParams(subjectID,trialID).torsoGamma = torsoGamma;
       subjectParams(subjectID,trialID).footm = footm;
       subjectParams(subjectID,trialID).footIxx = footIxx;
       subjectParams(subjectID,trialID).footIyy = footIyy; 
       subjectParams(subjectID,trialID).footIzz = footIzz;
       subjectParams(subjectID,trialID).torsom = torsom;
       subjectParams(subjectID,trialID).torsoIxx = torsoIxx;
       subjectParams(subjectID,trialID).torsoIyy = torsoIyy;
       subjectParams(subjectID,trialID).torsoIzz = torsoIzz;
       subjectParams(subjectID,trialID).legm = legm;
       subjectParams(subjectID,trialID).legIxx = legIxx;
       subjectParams(subjectID,trialID).legIyy = legIyy;
       subjectParams(subjectID,trialID).legIzz = legIzz;
       
       
    end
end

save('./experiments/humanFixedBase/subjectSizeParams.mat','subjectParams');




