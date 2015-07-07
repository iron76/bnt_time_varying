cl% load data

load('IMU_VICON_ShiftedData.mat');

isTest = 'false';

subjectList = 1:3;
trialList = 1:2 ; 

%subjectList = 3;
%trialList = 1 ; 

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d\nTrial : ',subjectID);
    for trialID = trialList
        fprintf('%d, ',trialID);
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
        
        P_PWA_C = temp.P_PWA_C(1:pSelec,:);
        fx_PWAPWA_1 = temp.fx_PWAPWA_1(1:pSelec,:);
        fx_PWAPWA_1(:,4:6) = fx_PWAPWA_1(:,4:6).*1e-3;
        
        P_G_1 = computeCentroidOfPoints(P_G_lankle,P_G_rankle);
        P_G_2 = computeCentroidOfPoints(P_G_lhip,P_G_rhip);
        P_G_3 = computeCentroidOfTriangle(P_G_lsho,P_G_rsho,P_G_tors);
        [R_G_0,P_G_0] = computeFootRotation(P_G_lhee,P_G_rhee,P_G_ltoe,P_G_rtoe); 
      
        R_0_G = R_G_0';
        R_G_1 = R_G_0;     
       
        R_1_G = R_G_1';
        
        len = size(P_G_1,1);
        %x_1 = repmat([1,0,0],len,1);
        x_1 = repmat([1,0,0],len,1);
        %% recomputing wrt to x-axis of frame1 to obtain both positive and negative angles
        
        %y_1_0 = repmat([0,-1,0],len,1);
        q1 = computeAngleBetweenVectors(x_1,(R_1_G*(P_G_2-P_G_1)')');
        q1 = repmat(0.5*pi,size(q1))-q1;
         %q1 = computeAngleBetweenVectors(y_1_0,(P_G_2-P_G_1));
        subplot(2,1,1);
        plot(temp.t_vicon(:,1:pSelec),q1.*(180/pi),'r'); hold on;
        xlabel('Time t(sec)');
        ylabel('q_1 and q_2 (degrees)');
        P_2_3 = recomputeVectorInFrame2(q1,R_1_G,(P_G_3-P_G_2));
        %q2 = computeAngleBetweenVectors((R_G_1*(P_G_2-P_G_1)')',(R_G_1*(P_G_3-P_G_2)')');
        x_2 = repmat([0,1,0],len,1);%x_1;
        q2 = computeAngleBetweenVectors(x_2,(P_G_3-P_G_2));
        q2 =  -repmat(0.5*pi,size(q2)) + q2;
        plot(temp.t_vicon(:,1:pSelec),q2.*(180/pi));
        legend('q_1 angle','q_2 angle');
        axis tight;
        
        subplot(2,1,2);
        q1s = sgolayfilt_wrapper(q1.*(180/pi),3,57);
        plot(temp.t_vicon(:,1:pSelec),q1s,'r'); hold on;
        xlabel('Time t(sec)');
        ylabel('q_1 and q_2 (degrees)');
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt_wrapper(q2.*(180/pi),3,57));
        legend('q_1 angle','q_2 angle');
        axis tight;
        
        q1 = sgolayfilt_wrapper(q1,3,57);
        q2 = sgolayfilt_wrapper(q2,3,57);
        
        dq1 = diff(q1)./1e-3;%diff(temp.t_vicon(1:end-1));
        dq2 = diff(q2)./1e-3;%diff(temp.t_vicon(1:end-1));
        dq1 = [dq1;dq1(end,:)];
        dq2 = [dq2;dq2(end,:)];
        
        dq1 = sgolayfilt_wrapper(dq1,3,57);
        dq2 = sgolayfilt_wrapper(dq2,3,57);
        
        ddq1 = diff(dq1)./1e-3;
        ddq2 = diff(dq2)./1e-3;
        
        ddq1 = sgolayfilt_wrapper(ddq1,3,57);
        ddq2 = sgolayfilt_wrapper(ddq2,3,57);
        
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
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt_wrapper(dq1.*(180/pi),3,57),'r'); hold on;
        xlabel('Time t(sec)');
        ylabel('dq_1 and dq_2 (degrees/sec)');
        
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt_wrapper(dq2.*(180/pi),3,57));
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
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt_wrapper(ddq1.*(180/pi),3,57),'r'); hold on;
        xlabel('Time t(sec)');
        ylabel('ddq_1 and ddq_2 (degrees/sec^2)');
        
        plot(temp.t_vicon(:,1:pSelec),sgolayfilt_wrapper(ddq2.*(180/pi),3,57));
        legend('ddq_1','ddq_2');
        axis tight;
       
        % computing R_G_imu
        [R_G_imuini,P_G_imuini] = computeInitialIMURotation(P_G_imuA,P_G_imuB,P_G_imuC);
        
        R_G_PWA = [1 0 0 ;0 -1 0;0 0 -1]; 
        %CLA: rotation matrix from PWA to G frame
        
        R_0_PWA = R_0_G * R_G_PWA; 
        %CLA: rotation matrix from PWA to 0 frame as the composition of two
        %R matrices
        
        r_G_from0toPWA = computeVectorFromPoints(repmat(-P_G_0,size(P_PWA_C,1),1),(R_G_PWA*P_PWA_C')')*1e-3;% positions in mm
        r_0_from0toPWA = (R_0_G * r_G_from0toPWA')'; 
        
        fx_0_1 = zeros(size(fx_PWAPWA_1));
        a_2_imulin = zeros(size(fx_PWAPWA_1,1),3);%
        v_2_imurot = zeros(size(fx_PWAPWA_1,1),3);%
        
        
       % R_G_imu0fake = [0 -1 0; 0 0 -1;1 0 0];
        R_0_1ini = euler2dcm([0,mean(q1(1:10)),0]);
                %R_0_1cla{i} = computeRy(q1(i));
         
        R_1_2ini = euler2dcm([0,mean(q2(1:10)),0]);
                %R_1_2cla{i} = computeRy(q2(i));

        R_G_2ini = R_G_0 * R_0_1ini * R_1_2ini;
       % R_G_imuinifake = [0 0 1; 0 -1 0; -1 0 0];
        R_2_imuini = R_G_2ini'* R_G_imuini;
        %R_2_imuini = R_G_2ini'* R_G_imuinifake;
        
       % R_G_0
        
       % R_0_1ini
       % R_1_2ini
        
        figure;
        subplot(2,1,1);
        plot(temp.t_imu,temp.a_imu_imulin');axis tight;
        xlabel('time (sec)');
        ylabel('m/sec^2');
        legend('x','y','z');
        title('Raw Acceleration of link2 (imuframe)');
        subplot(2,1,2);
        plot(temp.t_imu,temp.v_imu_imurot');  axis tight;
        xlabel('time (sec)');
        ylabel('rad/sec');
        legend('x','y','z');
        title('Raw Angular Velocity of link2 (imuframe)');
        
        
          for i = 1:length(temp.t_vicon)   
                a_2_imulin(i,:) =  (R_2_imuini*temp.a_imu_imulin(i,:)')'; %
                v_2_imurot(i,:) =  (R_2_imuini*temp.v_imu_imurot(i,:)')'; %

               % adjT_0_PWA{i} = [ R_0_PWA , zeros(3) ; -skew(r_0_from0toPWA(i,:)') * R_0_PWA , R_0_PWA]; 
                adjT_0_PWA{i} = [  R_0_PWA , -skew(r_0_from0toPWA(i,:)') * R_0_PWA ; zeros(3) , R_0_PWA  ]; 
                %CLA: the notation for force transformation is modified because
                %whe use the notation linear-angular in 6d vectors and not
                %angular-linear like in the Featherstone.

                fx_0_1(i,:) = ((adjT_0_PWA{i}) * fx_PWAPWA_1(i,:)')';
          end
        
%         R_0_1CELL = cell2mat(R_0_1);
%         R_0_1claCELL = cell2mat(R_0_1cla);
%         if (sum(R_0_1CELL-R_0_1claCELL) ~= 0);
%         disp('Something wrong!');
%         else disp('Two methods for R computation are equivalent.');
%         end
%         
%         R_1_2CELL = cell2mat(R_1_2);
%         R_1_2claCELL = cell2mat(R_1_2cla);
%         if (sum(R_1_2CELL-R_1_2claCELL) ~= 0);
%         disp('Something wrong!');
%         else disp('Two methods for R computation are equivalent.');
%         end
        
        figure;
        subplot(2,1,1);
        plot(temp.t_vicon,a_2_imulin); axis tight;
        xlabel('time (sec)');
        ylabel('m/sec^2');
        legend('x','y','z');
        title('Acceleration of link2');
        subplot(2,1,2);
        plot(temp.t_vicon,v_2_imurot);  axis tight;
        xlabel('time (sec)');
        ylabel('rad/sec');
        legend('x','y','z');
        title('Angular Velocity of link2');
        
        figure;
        subplot(2,1,1);
        plot(temp.t_vicon,fx_PWAPWA_1(:,1:3)); axis tight;
        xlabel('time (sec)');
        ylabel('Wrench Force(N)');
        title('Wrench measured');
        legend('x','y','z');
        subplot(2,1,2);
        plot(temp.t_vicon,fx_PWAPWA_1(:,4:6)); axis tight;
        xlabel('time (sec)');
        ylabel('Wrench Moments measured (Nm)');
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
        ylabel('Wrench Moments (Nm)');
        legend('x','y','z');
        
        processedSensorData(subjectID,trialID).R_G_0 = R_G_0;
        processedSensorData(subjectID,trialID).R_0_1 = R_0_1ini;
        processedSensorData(subjectID,trialID).R_2_imu = R_2_imuini;
        processedSensorData(subjectID,trialID).R_G_imu0 = R_G_imuini;
        processedSensorData(subjectID,trialID).R_G_2 = R_G_2ini;
        processedSensorData(subjectID,trialID).R_G_PWA = R_G_PWA;
        processedSensorData(subjectID,trialID).R_0_PWA = R_0_PWA;
        processedSensorData(subjectID,trialID).q1 = q1;
        processedSensorData(subjectID,trialID).q2 = q2;
        processedSensorData(subjectID,trialID).dq1 = dq1;
        processedSensorData(subjectID,trialID).dq2 = dq2;
        processedSensorData(subjectID,trialID).ddq1 = ddq1;
        processedSensorData(subjectID,trialID).ddq2 = ddq2;
        processedSensorData(subjectID,trialID).a_2_imulin = a_2_imulin;
        processedSensorData(subjectID,trialID).v_2_imurot = v_2_imurot;
        processedSensorData(subjectID,trialID).adjT_0_PWA = adjT_0_PWA;
        processedSensorData(subjectID,trialID).f_x_1 = fx_0_1';
%         processedSensorData(subjectID,trialID).a_2_imulin = a_2_imulin;
%         processedSensorData(subjectID,trialID).v_2_imurot = v_2_imurot;
        processedSensorData(subjectID,trialID).t = temp.t_vicon;
        processedSensorData(subjectID,trialID).imu = [a_2_imulin v_2_imurot]';
        %processedSensorData(subjectID,trialID).imu = 
        processedSensorData(subjectID,trialID).ftx = fx_0_1';
        
    end
end
% make loop
%processedSensorData = imu_vicon_shiftedData;
if(strcmp(isTest,'true')~=1)
    save('./experiments/humanFixedBase/preProcessedSensorData.mat','processedSensorData');%a_2_imulin
end

