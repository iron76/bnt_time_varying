
% ======loading data from Vicon
load('IMU_VICON_ShiftedData.mat');

isTest = 'false';

subjectList = 1;
trialList = 1 ; 

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d\nTrial : ',subjectID);
    for trialID = trialList
        fprintf('%d, ',trialID);
        
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
        
        f_fp = temp.f_fp(1:pSelec,:);
        f_fp(:,1:3) = temp.f_fp(:,1:3)*1e-3; % converting moments from Nmm to Nm
 
        
        P_G_1 = computeCentroidOfPoints(P_G_lankle,P_G_rankle);
        P_G_2 = computeCentroidOfPoints(P_G_lhip,P_G_rhip);
        P_G_3 = computeCentroidOfTriangle(P_G_lsho,P_G_rsho,P_G_tors);
        
        [R_0_G,P_G_0] = computeFootRotation(P_G_lhee,P_G_rhee,P_G_ltoe,P_G_rtoe); 
        
        R_G_0 = R_0_G';
        R_G_1 = R_G_0;     % because point P_G_1 is fixed on the foot in URDF
        

        
        %% Computing q1 and q2  angles 
       
        len = size(P_G_1,1);
        
        % ======JOINT ANGLE q1
        l1 = (P_G_2 - P_G_1);
        q1 = zeros (len, 1);
        
        for i = 1 : len;
            q1(i) =atan2(-l1(i,2),l1(i,3));
        end
       
        % ======JOINT ANGLE q2
        l2 = (P_G_3-P_G_2);
        q_temp = zeros (len, 1);

        for i = 1 : len;
            q_temp(i) = atan2(-l2(i,2),l2(i,3));
        end
        
        q2 = q_temp-q1;
    
        
        %% Using Savitzky-Golay filtering for differentiation

        window = 57;
        [~, diffCoeff] = sgolay_wrapper(3, window);
        %diffCoeff is a matrix of (polynomialOrder-1) columns where:
        %- ( ,1) --> coefficient for S-Golay as smoother;
        %- ( ,2) --> coefficient for S-Golay as 1st differentiator;
        %- ( ,3) --> coefficient for S-Golay as 2nd differentiator;
        %  .   
        %  .
        %  .
        %- ( ,polynomialOrder-1) --> coefficient for S-Golay as (polynomialOrder) differentiator;
        
        halfWindow  = ((window+1)/2) -1;
        dq1_sg = zeros(len, 1);
        dq2_sg = zeros(len, 1);
        ddq1_sg = zeros(len, 1);
        ddq2_sg = zeros(len, 1);
        
        for n = (window+1)/2:len-(window+1)/2,
              % 1st differential
              dq1_sg(n) = dot(diffCoeff(:,2),q1(n - halfWindow:n + halfWindow));
              dq2_sg(n) = dot(diffCoeff(:,2),q2(n - halfWindow:n + halfWindow));
              % 2nd differential
              ddq1_sg(n) = dot(diffCoeff(:,3),q1(n - halfWindow:n + halfWindow));
              ddq2_sg(n) = dot(diffCoeff(:,3),q2(n - halfWindow:n + halfWindow));
        end
        
        deltaT = 1e-3;
        dq1_sg = dq1_sg ./deltaT;
        dq2_sg = dq2_sg ./ deltaT;
        ddq1_sg = ddq1_sg ./ (deltaT)^2;
        ddq2_sg = ddq2_sg ./ (deltaT)^2;
        
       
        figure;
        subplot(311);
        plot1 = plot(temp.t_vicon(:,1:pSelec),q1.*(180/pi),'lineWidth',1.0); hold on;
        set(plot1,'color',[1 0 0]);
        plot2= plot(temp.t_vicon(:,1:pSelec),q2.*(180/pi),'lineWidth',1.0); hold on;
        set(plot2,'color',[0 0.498039215803146 0]);
        leg = legend('$q_1$','$q_2$','Location','northeast');
        title('Joint Quantities','FontSize',15);
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',15);
        xlabel('Time [s]','FontSize',15);
        ylabel('Angle [deg]','FontSize',15);
        axis tight;
        grid on;  

        subplot(312);
        plot1 = plot(temp.t_vicon(:,1:pSelec),(180/pi)*dq1_sg,'lineWidth',1.0); hold on;
        set(plot1,'color',[1 0 0]);
        plot2= plot(temp.t_vicon(:,1:pSelec),(180/pi)*dq2_sg,'lineWidth',1.0); hold on;
        set(plot2,'color',[0 0.498039215803146 0]);
        leg = legend('$\dot q_{1}$','$\dot q_{2}$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',15);
        xlabel('Time [s]','FontSize',15);
        ylabel('Velocity [deg/s]','FontSize',15);
        axis tight;
        grid on;
        
        subplot(313);
        plot1 = plot(temp.t_vicon(:,1:pSelec),(180/pi)*ddq1_sg,'lineWidth',1.0); hold on;
        set(plot1,'color',[1 0 0]);
        plot2= plot(temp.t_vicon(:,1:pSelec),(180/pi)*ddq2_sg,'lineWidth',1.0); hold on;
        set(plot2,'color',[0 0.498039215803146 0]);
        leg = legend('$\ddot q_{1}$','$\ddot q_{2}$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',15);
        xlabel('Time [s]','FontSize',15);
        ylabel('Acceleration [deg/s^2]','FontSize',15);
        axis tight;
        grid on;
        
       
        %% IMU sensor
      
        % =======Plotting raw data coming from IMU sensor in IMU frame
        
        figure;
        subplot(211);
        plot(temp.t_imu,temp.a_imu_imulin'); axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Linear Acceleration [m/sec^2]','FontSize',15);
        legend('$a_x$','$a_y$','$a_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',15);
        title('Raw IMU data of link 2 (IMU frame)','FontSize',15);
        grid on;
        
        subplot(212);
        plot(temp.t_imu,temp.v_imu_imurot');  axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Angular Velocity [rad/s]','FontSize',15);
        legend('$w_x$','$w_y$','$w_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',15);
        grid on;
        
        
        % ======Extracting rotation matrix R_2_imuini
        
       R_0_1ini = euler2dcm([0,mean(q1(1:10)),0]);
       R_1_2ini = euler2dcm([0,mean(q2(1:10)),0]); 
       R_G_2ini = R_G_0 * R_0_1ini * R_1_2ini;
        
       [R_G_imuini,~] = computeInitialIMURotation(P_G_imuA,P_G_imuB,P_G_imuC);
       R_2_imuini = R_G_2ini'* R_G_imuini;
%   
%         R_1ini_0 = euler2dcm([0,mean(q1(1:10)),0]); 
%         R_2ini_1 = euler2dcm([0,mean(q2(1:10)),0]); 
%         R_2ini_G =R_2ini_1 * R_1ini_0 * R_0_G
%         
%         [R_G_imuini,~] = computeInitialIMURotation(P_G_imuA,P_G_imuB,P_G_imuC);
%         R_2_imuini = R_imuini_G' * R_2ini_G';
  
     
        
        % ======Computing IMU data in link 2 frame
        
        a_2_imulin = zeros(size(q1,1),3);
        v_2_imurot = zeros(size(q1,1),3);
        
          for i = 1:length(temp.t_vicon)   
                a_2_imulin(i,:) =  (R_2_imuini*temp.a_imu_imulin(i,:)')'; 
                v_2_imurot(i,:) =  (R_2_imuini*temp.v_imu_imurot(i,:)')'; 
          end
        
          
        % ======Plotting data coming from IMU sensor in link 2 frame
          
        figure;
        subplot(211);
        plot(temp.t_vicon,a_2_imulin); axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Linear Acceleration [m/sec^2]','FontSize',15);
        legend('$a_x$','$a_y$','$a_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',15);
        title('IMU data of link 2 (link 2 frame)','FontSize',15);
        grid on;
        
        subplot(212);
        plot(temp.t_vicon,v_2_imurot);  axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Angular Velocity [rad/s]','FontSize',15);
        legend('$w_x$','$w_y$','$w_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',15);
        grid on;
            
        
        %% Force plate sensing
        
        %======fixed rotation from Global and Force plate reference frames
        R_G_fp = [-1 0  0;
                   0 1  0; 
                   0 0 -1];   
              
        R_0_fp = R_0_G * R_G_fp; 
        
        % ======center of force plate in mm (below the force plate) in Global frame
        P_G_fp = [231.75,254,-43.3]; 
        
        r_G_from0toFp = P_G_0 - P_G_fp;
        r_0_from0toFp = R_0_G*r_G_from0toFp';
        
        r_0_from0toFpm = r_0_from0toFp*1e-3; %converting to m

        f_0 = zeros(length(temp.t_vicon),6);
        XStar_0_fp = [R_0_fp' skew(r_0_from0toFpm)*R_0_fp'; zeros(3) R_0_fp'];
        f_0 = (XStar_0_fp * f_fp')';
        
       
        figure;
        subplot(211);
        plot(temp.t_vicon,f_fp(:,4:6)); axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Force [N]','Fontsize',15);
        title('Wrench measured in force plate (Fp) frame','FontSize',15);
        legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
        grid on;
        
        subplot(212);
        plot(temp.t_vicon,f_fp(:,1:3)); axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Moment [Nm]','FontSize',15);
        legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
        grid on;
        
       
        figure;
        subplot(211);
        plot(temp.t_vicon, f_0(:,4:6)); axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Force [N]','FontSize',15);
        title('Wrench measured in link 0 frame','FontSize',15);
        legend('$F_x$','$F_y$','$F_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
        grid on;
        
        subplot(212);
        plot(temp.t_vicon, f_0(:,1:3)); axis tight;
        xlabel('Time [s]','FontSize',15);
        ylabel('Moment [Nm]','FontSize',15);
        legend('$M_x$','$M_y$','$M_z$','Location','northeast');
        set(legend,'Interpreter','latex');
        set(legend,'FontSize',20);
        grid on;
        
        
        
        %% clustering data
        
        processedSensorData(subjectID,trialID).R_G_0 = R_G_0;
        processedSensorData(subjectID,trialID).R_0_1 = R_0_1ini; 
        processedSensorData(subjectID,trialID).R_2_imu = R_2_imuini; 
        processedSensorData(subjectID,trialID).R_G_imu0 = R_G_imuini; 
        processedSensorData(subjectID,trialID).R_G_2 = R_G_2ini; 
        processedSensorData(subjectID,trialID).R_G_fp = R_G_fp;
        processedSensorData(subjectID,trialID).R_0_fp = R_0_fp;
        processedSensorData(subjectID,trialID).q1 = q1;
        processedSensorData(subjectID,trialID).q2 = q2;
        processedSensorData(subjectID,trialID).dq1 = dq1_sg;
        processedSensorData(subjectID,trialID).dq2 = dq2_sg;
        processedSensorData(subjectID,trialID).ddq1 = ddq1_sg;
        processedSensorData(subjectID,trialID).ddq2 = ddq2_sg;
        processedSensorData(subjectID,trialID).a_2_imulin = a_2_imulin';
        processedSensorData(subjectID,trialID).v_2_imurot = v_2_imurot';
        processedSensorData(subjectID,trialID).XStar_0_fp = XStar_0_fp;
        processedSensorData(subjectID,trialID).t = temp.t_vicon;
        processedSensorData(subjectID,trialID).imu = [a_2_imulin v_2_imurot]';
        processedSensorData(subjectID,trialID).f_0 = f_0';
       
        
    end
end

if(strcmp(isTest,'true')~=1)
    save('./experiments/humanFixedBase/data/preProcessedSensorData.mat','processedSensorData');
end

