% load data

load('IMU_VICON_ShiftedData.mat');

%[footAlpha.footBeta,footGamma]

subjectList = 1;
trialList = 1 ; 

for subjectID = 1:length(subjectList)
    for trialID = 1:length(trialList)
        figure;
        temp = imu_vicon_shiftedData(subjectID,trialID);
        
        
        pSelec = size(temp.P_G_lhee,1);
        
        P_G_lhee = temp.P_G_lhee(1:pSelec,:);
        P_G_ltoe = temp.P_G_ltoe(1:pSelec,:);
        P_G_rhee = temp.P_G_rhee(1:pSelec,:);
        P_G_rtoe = temp.P_G_rtoe(1:pSelec,:);
        
        P_G_lankle = temp.P_G_lankle(1:pSelec,:);
        P_G_lhip = temp.P_G_lhip(1:pSelec,:);
        P_G_lsho = temp.P_G_lsho(1:pSelec,:);
        
        P_G_rankle = temp.P_G_rankle(1:pSelec,:);
        P_G_rhip = temp.P_G_rhip(1:pSelec,:);
        P_G_rsho = temp.P_G_rsho(1:pSelec,:);
        
        P_G_tors = temp.P_G_tors(1:pSelec,:);
        
        P_G_imuA = temp.P_G_imuA(1:pSelec,:);
        P_G_imuB = temp.P_G_imuB(1:pSelec,:);
        P_G_imuC = temp.P_G_imuC(1:pSelec,:);
        
        P_G_s = temp.P_G_s(1:pSelec,:);
        f_GPWA_s = temp.f_GPWA_s(1:pSelec,:);
        f_GPWA_s(:,4:6) = f_GPWA_s(:,4:6).*1e-3;
        
        P_G_1 = computeCentroidOfPoints(P_G_lankle,P_G_rankle);
        P_G_2 = computeCentroidOfPoints(P_G_lhip,P_G_rhip);
        P_G_3 = computeCentroidOfTriangle(P_G_lsho,P_G_rsho,P_G_tors);
        
        [R_G_0, P_G_0] = computeFootRotation(P_G_lhee,P_G_rhee,P_G_ltoe,P_G_rtoe);    
        R_0_G = R_G_0;
        R_G_1 = R_G_0;        
        P_0_0PWA = computeVectorFromPoints(repmat(P_G_0,size(P_G_s,1),1),P_G_s).*1e-3;
        len = size(P_G_1,1);
        %z_1_0 = repmat((R_G_0*[0,0,1]')',len,1);
        z_1_0 = repmat([0,0,1],len,1);
         %q1 = computeAngleBetweenVectors(z_1_0,(R_G_1*(P_G_2-P_G_1)')');
        q1 = computeAngleBetweenVectors(z_1_0,P_G_2-P_G_1);
        subplot(2,1,1);
        plot(temp.t_vicon(:,1:pSelec),q1.*(180/pi),'r'); hold on;
        xlabel('Time t(sec)');
        ylabel('q_1 and q_2 (degrees)');
        %q2 = computeAngleBetweenVectors((R_G_1*(P_G_2-P_G_1)')',(R_G_1*(P_G_3-P_G_2)')');
        q2 = computeAngleBetweenVectors((R_G_1*(P_G_3-P_G_2)')',(R_G_1*(P_G_2-P_G_1)')');
        plot(temp.t_vicon(:,1:pSelec),q2.*(180/pi));
        legend('q_1','q_2');
        axis tight;
        
        subplot(2,1,2);
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt(q1.*(180/pi),3,57),'r'); hold on;
        xlabel('Time t(sec)');
        ylabel('q_1 and q_2 (degrees)');
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt(q2.*(180/pi),3,57));
        legend('q_1','q_2');
        axis tight;
        
        q1 = sgolayfilt(q1,3,57);
        q2 = sgolayfilt(q2,3,57);
        
        dq1 = diff(q1)./1e-3;%diff(temp.t_vicon(1:end-1));
        dq2 = diff(q2)./1e-3;%diff(temp.t_vicon(1:end-1));
        dq1 = [dq1;dq1(end,:)];
        dq2 = [dq2;dq2(end,:)];
        
        dq1 = sgolayfilt(dq1,3,57);
        dq2 = sgolayfilt(dq2,3,57);
        
        ddq1 = diff(dq1)./1e-3;
        ddq2 = diff(dq2)./1e-3;
        
        ddq1 = sgolayfilt(ddq1,3,57);
        ddq2 = sgolayfilt(ddq2,3,57);
        
        ddq1 = [ddq1;ddq1(end,:)];
        ddq2 = [ddq2;ddq2(end,:)];
        
        
        figure;
        subplot(2,1,1);
        plot(temp.t_vicon(:,1:pSelec),dq1.*(180/pi),'r'); hold on;
        xlabel('Time t(sec)');
        ylabel('dq_1 and dq_2 (degrees/sec)');
        
        plot(temp.t_vicon(:,1:pSelec),dq2.*(180/pi));
        legend('dq_1','dq_2');
        axis tight;
        
        subplot(2,1,2);
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt(dq1.*(180/pi),3,57),'r'); hold on;
        xlabel('Time t(sec)');
        ylabel('dq_1 and dq_2 (degrees/sec)');
        
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt(dq2.*(180/pi),3,57));
        legend('dq_1','dq_2');
        axis tight;
        
        figure;
        subplot(2,1,1);
        plot(temp.t_vicon(:,1:pSelec),ddq1.*(180/pi),'r'); hold on;
        xlabel('Time t(sec)');
        ylabel('ddq_1 and ddq_2 (degrees/sec^2)');
        
        plot(temp.t_vicon(:,1:pSelec),ddq2.*(180/pi));
        legend('ddq_1','ddq_2');
        axis tight;
        subplot(2,1,2);
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt(ddq1.*(180/pi),3,57),'r'); hold on;
        xlabel('Time t(sec)');
        ylabel('ddq_1 and ddq_2 (degrees/sec^2)');
        
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt(ddq2.*(180/pi),3,57));
        legend('ddq_1','ddq_2');
        axis tight;
        % computing R_G_imu
        
        [R_G_imu0,p_G_imu0] = computeInitialIMURotation(P_G_imuA,P_G_imuB,P_G_imuC);
        
        
        fx_0_1 = zeros(size(f_GPWA_s));
        a_2_imulin = zeros(size(f_GPWA_s,1),3);
        v_2_imurot = zeros(size(f_GPWA_s,1),3);
        
        for i = 1:length(temp.t_vicon)
            R_0_1{i} = euler2dcm([0,0,q1(i)]);
            R_1_2{i} = euler2dcm([0,0,q2(i)]);
            
            R_G_2{i} = R_G_0 * R_0_1{i}  * R_1_2{i};
            R_2_imu{i} = R_G_2{i}'*R_G_imu0;
            
            a_2_imulin(i,:) =  (R_2_imu{i}*temp.a_imu_imulin(i,:)')';
            v_2_imurot(i,:) =  (R_2_imu{i}*temp.v_imu_imurot(i,:)')'; 
            adjT_0_PWAG{i} = [ R_0_G , zeros(3) ; -skew(P_0_0PWA(i,:)') * R_0_G , R_0_G]; 
            
            fx_0_1(i,:) = (adjT_0_PWAG{i}'*f_GPWA_s(i,:)')';
        end
        
        figure;
        subplot(2,1,1);
        plot(temp.t_vicon,a_2_imulin); axis tight;
        xlabel('time (sec)');
        ylabel('m/sec^2');
        legend('x','y','z');
        title('Acceleration of link2');
        subplot(2,1,2);
        plot(temp.t_vicon,v_2_imurot);  axis tight;
        title('AngularVelocity of link2');
        xlabel('time (sec)');
        ylabel('rad/sec');
        legend('x','y','z');
        
        figure;
        subplot(2,1,1);
        plot(temp.t_vicon,f_GPWA_s(:,1:3)); axis tight;
        xlabel('time (sec)');
        ylabel('Wrench Force(N)');
        title('Wrench measured');
        legend('x','y','z');
        subplot(2,1,2);
        plot(temp.t_vicon,f_GPWA_s(:,4:6)); axis tight;
        xlabel('time (sec)');
        ylabel('Wrench Momments(Nm)');
        legend('x','y','z');
        
        figure;
        subplot(2,1,1);
        plot(temp.t_vicon, fx_0_1(:,1:3)); axis tight;
        xlabel('time (sec)');
        ylabel('Wrench Force (N)');
        title('External Wrench at link0');
        legend('x','y','z');
        subplot(2,1,2);
        plot(temp.t_vicon, fx_0_1(:,4:6)); axis tight;
        xlabel('time (sec)');
        ylabel('Wrench Momments (Nm)');
        legend('x','y','z');
        
        processedSensorData(subjectID,trialID).R_G_0=R_G_0;
        processedSensorData(subjectID,trialID).R_2_imu = R_2_imu;
        processedSensorData(subjectID,trialID).R_G_imu0 = R_G_imu0;
        processedSensorData(subjectID,trialID).q1 = q1;
        processedSensorData(subjectID,trialID).q2 = q2;
        processedSensorData(subjectID,trialID).dq1 = dq1;
        processedSensorData(subjectID,trialID).dq2 = dq2;
        processedSensorData(subjectID,trialID).ddq1 = ddq1;
        processedSensorData(subjectID,trialID).ddq2 = ddq2;
        
        processedSensorData(subjectID,trialID).a_2_imulin = a_2_imulin;
        processedSensorData(subjectID,trialID).v_2_imurot = v_2_imurot;
        processedSensorData(subjectID,trialID).adjT_0_PWA = adjT_0_PWAG;
        processedSensorData(subjectID,trialID).f_x_1 = fx_0_1';

        processedSensorData(subjectID,trialID).a_2_imulin = a_2_imulin;
        processedSensorData(subjectID,trialID).v_2_imurot = v_2_imurot;
        
        processedSensorData(subjectID,trialID).t = temp.t_vicon;
        processedSensorData(subjectID,trialID).imu = [a_2_imulin v_2_imurot]';
        %processedSensorData(subjectID,trialID).imu = 
        processedSensorData(subjectID,trialID).ftx = fx_0_1';
        
    end
end
% make loop
%processedSensorData = imu_vicon_shiftedData;
save('./experiments/humanFixedBase/preProcessedSensorData.mat','processedSensorData');%a_2_imulin