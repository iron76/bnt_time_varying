% predictSensorMeasFromRNEA
% Script makes a prediction of the sensor measurements using variables from
% RNEA computation. 
%
% Toolbox requirements:
% - iDynTree - mex
% - Featherstone toolbox (v2)

% Measurements equations are:
% y_2(acc)  = S_lin*(imu_X_2 * a_2_2) + S_ang*(imu_X_2 * v_2_2) x S_lin*(imu_X_2 * v_2_2)
% y_2(gyro) = S_ang*(imu_X_2 * v_2_2)
% y_1(fp)   = fp_XStar_0 * ((0_XStar_1*f_1_1)) - I_0 * ag)
%
% To predict  measurements we use quantities coming from RNEA computation.


clear;clc;close all

%% testOptions
plotJointQuantities = false;
plotLinkQuantities = false; 
plotSensorPrediction = true; 
plotResultsVariances = true;
%analyseRNEA = false; % To run checkRNEA against iDynTree.Not yet tested with new refactor (pls retain to false)

%% selected subjects and trials
subjectList = 1;
trialList = 1;  

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d',subjectID);
    for trialID = trialList
         fprintf('\nTrial : %d\n',trialID);
         
        %% load data 
        load('./experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat');

        currentTrial = synchronisedData(subjectID,trialID);

        dataTime = currentTrial.dataTime;
        len = length(dataTime);
        
        q = currentTrial.q;
        dq = currentTrial.dq;
        ddq = currentTrial.ddq;

        wrench_fp_fp = currentTrial.wrench_fp_fp;
        imu = currentTrial.imu;
        
        %% build models

        load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');
        
        
        % =====structure of sensors
        sens.parts    = {'leg'         ,'torso'};       %force of the forceplate is trasmitted into the leg
        sens.labels   = {'fts'         ,'imu'  };
        sens.ndof     = {6             ,6      };
        
        currentModel = humanThreeLinkModelFromURDF(subjectID);
        dmodel  = currentModel.dmodel;                      %deterministic model 
        ymodel  = humanThreeLinkSens(dmodel, sens);  

        dmodel  = autoTreeStochastic(dmodel, 1e-6);         % probabilistic model for D equation (added Sv and Sw)
        ymodel  = humanThreeLinkSensStochastic(ymodel);     % probabilistic model for Y(q,dq) d = y (added Sy)

        myModel = model(dmodel);
        mySens  = sensors(ymodel);  
        
        %% RNEA class for computing prediction variables
        
        ymodel_RNEA  = autoSensRNEA(dmodel);
        mySens_RNEA  = sensors(ymodel_RNEA);
        myRNEA       = RNEA(myModel, mySens_RNEA);

        y_RNEA_f = zeros(6*dmodel.NB, len);
        y_RNEA_ddq = zeros(dmodel.NB, len);
        fx = cell(dmodel.NB);

        % Ordering y_RNEA in the form [fx1 fx2 ddq1 ddq2]
        for i = 1 : dmodel.NB
            for t = 1 : len
                fx{i,1} = zeros(6,1); 
                y_RNEA_f(((1:6)+(i-1)*6), t) = [fx{i,1}];
                y_RNEA_ddq(i, t) = [ddq(t,i)];
            end
            y_RNEA = [y_RNEA_f ; y_RNEA_ddq];
        end

        
        v_2_RNEA = zeros(len,6);
        a_2_RNEA = zeros(len,6);
        f_1_RNEA = zeros(len,6);
        tau_RNEA = zeros(length(dataTime),dmodel.NB);
        
        for i = 1 : len
             myRNEA = myRNEA.setState(q(i,:)', dq(i,:)');
             myRNEA = myRNEA.setY(y_RNEA(:,i));
             myRNEA = myRNEA.solveID();
             
%              v_RNEA{i,:} = myRNEA.v;      %containing both v_1_1 and v_2_2
%              a_RNEA{i,:} = myRNEA.a;      %containing both a_1_1 and a_2_2
%              f_RNEA{i,:} = myRNEA.f;      %containing both f_1_1 and f_2_2
%                        
             v_2_RNEA(i,:) = myRNEA.v(:,2);      %v_2_2
             a_2_RNEA(i,:) = myRNEA.a(:,2);      %a_2_2
             f_1_RNEA(i,:) = myRNEA.f(:,1);      %f_1_1 
             tau_RNEA(i,:) = myRNEA.tau; 
        end
        
        clear fx;
        % ========end RNEA
       
        %% Plot raw components
        
        if(plotJointQuantities)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;

            subplot(311);
            plot1 = plot(dataTime,(180/pi)*q(:,1),'lineWidth',2.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(dataTime,(180/pi)*q(:,2),'lineWidth',2.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$q_1$','$q_2$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Angle [deg]','FontSize',18);
            axis tight;
            grid on;

            subplot(312);
            plot1 = plot(dataTime,(180/pi)*dq(:,1),'lineWidth',2.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2= plot(dataTime,(180/pi)*dq(:,2),'lineWidth',2.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$\dot q_1$','$\dot q_2$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Velocity [deg/s]','FontSize',18);
            axis tight;
            grid on;

            subplot(313);
            plot1 = plot(dataTime,tau_RNEA(:,1), 'lineWidth',2.0); hold on;
            set(plot1,'color',[1 0 0]);
            plot2 = plot(dataTime,tau_RNEA(:,2), 'lineWidth',2.0); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);
            leg = legend('$\tau_1$','$\tau_2$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Torque [Nm]','FontSize',18);
            axis tight;
            grid on;

            axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
            text(0.5, 0.99,'\bf Joint Quantities','HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',14);
        end

        %%  ========================= Sensor prediction ===========================
        
        load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');

        currentTrialSens = sensorLinkTransforms(subjectID,trialID);

        X_imu_2 = currentTrialSens.X_imu_2;
        XStar_fp_0 = currentTrialSens.XStar_fp_0;
        XStar_0_1 = currentTrialSens.XStar_0_1;


        load('./experiments/humanFixedBase/data/subjectSizeParams.mat');

        currentParams = subjectParams(subjectID);

        footMass =  currentParams.footMass;
        posP_0 = [0; 0; (0.5*currentParams.footHeight)];
        footIxx =  currentParams.footIxx;
        footIyy =  currentParams.footIyy;
        footIzz =  currentParams.footIzz;

        S_lin = [zeros(3) eye(3)];
        S_ang = [eye(3) zeros(3)];


        fprintf('Predicting sensor measurement using RNEA computation\n');

        % ====IMU PREDICTION
        a_imu_imuPred = S_lin*(X_imu_2*a_2_RNEA') + skew(S_ang*X_imu_2*v_2_RNEA')*(S_lin*X_imu_2*v_2_RNEA'); %skew or cross product are equivalent
        omega_imu_imuPred = S_ang*(X_imu_2*v_2_RNEA');
        omegaDot_imu_imuPred = S_ang*(X_imu_2*a_2_RNEA');

        % ====FORCE PLATE PREDICTION
        a_G_grav = [0;0;0;0;0;-9.8100];    %Featherstone-like notation, in global reference
        X_G_0 = currentTrialSens.X_G_0;
        X_0_G = InverseAdjTransform(X_G_0);
        a_0_grav = X_0_G * a_G_grav;

        I_0 = createSpatialInertia(footIxx,footIyy,footIzz,footMass,posP_0);

        f_fp_fpPred = zeros(6,len);
        f_0_1 = zeros(6,len);
        f_1_RNEA = f_1_RNEA';

        for i = 1 : len
             f_0_1(:,i) = (XStar_0_1{i,1} * f_1_RNEA(:,i)); 
             f_fp_fpPred(:,i) = XStar_fp_0 * ((f_0_1(:,i))-(I_0 * a_0_grav));
        end
  
        %% Plot link quantities from ID Featherstone
        
        if(plotLinkQuantities)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;

            subplot(211);
            plot1 = plot(dataTime,S_lin*a_2_RNEA','lineWidth',2.0); hold on; 
            leg = legend('$a_x$','$a_y$','$a_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Lin Acceleration [m/sec^2]','FontSize',15);
            title('Lin Acceleration and Ang Velocity in frame associated to link2','FontSize',14);
            axis tight; 
            grid on;

            subplot(212);
            plot2 = plot(dataTime,S_ang*v_2_RNEA','lineWidth',2.0); hold on;
            leg = legend('$w_x$','$w_y$','$w_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Ang Velocity [rad/sec]','FontSize',15);
            axis tight; 
            grid on;

            figure()
            subplot(211)
            plot3 = plot(dataTime,f_1_RNEA(4:6,:),'lineWidth',2.0); hold on;
            leg = legend('$F_x$','$F_y$','$F_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Force [N]','FontSize',15);
            title('Wrench in frame associated to link1 ','FontSize',14);
            axis tight; 
            grid on;

            subplot(212)
            plot4 = plot(dataTime,f_1_RNEA(1:3,:),'lineWidth',2.0); hold on;
            leg = legend('$\mu_x$','$\mu_y$','$\mu_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Moment [Nm]','FontSize',15);
            axis tight; 
            grid on;    
        end

    %     if(analyseRNEA == true)
    %         checkRNEA_iDynTree;
    %     end
        %% Plot predictions 

        if(plotSensorPrediction)
        %% Accelerometer prediction comparison  --> linear acceleration

            fig = figure('name','RNEA');
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;

            subplot(221)
            plot1 = plot(dataTime,a_imu_imuPred(1,:),dataTime,a_imu_imuPred(2,:),dataTime,a_imu_imuPred(3,:), 'lineWidth',2.0);
            leg = legend('$a_x$','$a_y$','$a_z$','Location','southeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Prediction ','FontSize',15);
            title('Linear acceleration [m/sec^2]','FontSize',14);
            axis tight;
            grid on; 

            subplot(223)
            plot2 = plot(dataTime,imu(:,1),dataTime,imu(:,2),dataTime,imu(:,3),'lineWidth',2.0);
            leg = legend('$a_x$','$a_y$','$a_z$','Location','southeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Actual','FontSize',15);
            axis tight;
            grid on; 
        
            %% Gyroscope prediction comparison  --> angular velocities

            subplot(222)
            plot(dataTime,omega_imu_imuPred(1,:),dataTime,omega_imu_imuPred(2,:),dataTime,omega_imu_imuPred(3,:), 'lineWidth',2.0);
            leg = legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','southeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Prediction','FontSize',15);
            title('Angular velocity [rad/sec]','FontSize',14);
            axis tight;
            grid on; 

            subplot(224)
            plot(dataTime,imu(:,4),dataTime,imu(:,5),dataTime,imu(:,6),'lineWidth',2.0);
            leg = legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','southeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Actual','FontSize',15);
            axis tight;
            grid on; 

            axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
            text(0.5, 0.99,(sprintf('RNEA measurements prediction, sensor frame (Subject: %d, Trial: %d)',subjectID, trialID)),'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',16);

            %% Force prediction comparison  --> force

            fig = figure('name','RNEA');
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;

            subplot(221);
            plot(dataTime,f_fp_fpPred(4,:),dataTime,f_fp_fpPred(5,:),dataTime,f_fp_fpPred(6,:), 'lineWidth',2.0);
            leg = legend('$F_x$','$F_y$','$F_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Prediction','FontSize',15);
            title('Force [N]','FontSize',15);
            hold on; 
            axis tight;
            grid on;

            subplot(223)
            plot(dataTime,wrench_fp_fp(:,4),dataTime,wrench_fp_fp(:,5),dataTime,wrench_fp_fp(:,6),'lineWidth',2.0);
            leg = legend('$F_x$','$F_y$','$F_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Actual','FontSize',15);
            axis tight;
            grid on;

            %% Moment prediction comparison  --> moment 

            subplot(222);
            plot(dataTime,f_fp_fpPred(1,:),dataTime,f_fp_fpPred(2,:),dataTime,f_fp_fpPred(3,:), 'lineWidth',2.0);
            leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Prediction','FontSize',15);
            title('Moment [Nm]','FontSize',15);
            hold on;
            axis tight;
            grid on;

            subplot(224)
            plot(dataTime,wrench_fp_fp(:,1),dataTime,wrench_fp_fp(:,2),dataTime,wrench_fp_fp(:,3),'lineWidth',2.0);
            leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',15);
            ylabel('Actual','FontSize',15);
            axis tight;
            grid on;

            axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
            text(0.5, 0.99,(sprintf('RNEA  measurements prediction, sensor frame (Subject: %d, Trial: %d)',subjectID, trialID)),'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',14);

        end
    
        %% Plot variances
        
        if(plotResultsVariances)

            dataRNEA = [f_fp_fpPred; omegaDot_imu_imuPred; a_imu_imuPred];

            load('./experiments/humanFixedBase/intermediateDataFiles/finalResults.mat');
            currentMAP = finalResults(subjectID,trialID);

            data.Sy = currentMAP.resMAP.data.Sy;
            data.y =currentMAP.resMAP.data.y;

            for i = 1:12
                fig = figure();
                axes1 = axes('Parent',fig,'FontSize',16);
                box(axes1,'on');
                hold(axes1,'on');
                grid on;

                lineProps = {'LineWidth', 3.0};
                shad1 = shadedErrorBar(dataTime,data.y(i,:),sqrt(data.Sy(i,:)),lineProps);
                plot1 = plot(dataTime,dataRNEA(i,:),'lineWidth',1.0);hold on;
                set(plot1,'color',[1 0 0]);
                xlabel('Time [s]','FontSize',15);
                leg = legend([shad1.mainLine,shad1.patch, plot1], {'sensorData','sensorDataVariance','RNEAprediction'},'Location','southeast');
                title('Comparison between real data and RNEA prediction data','FontSize',13);

                if (i >=1 && i<4) 
                    ylabel('f1: Moment [Nm]','FontSize',15);
                elseif (i >=4 && i<7) 
                    ylabel('f1: Force[N]','FontSize',15);
                elseif (i >=7 && i<10)
                    ylabel('a2: Angular acceleration [rad/s]','FontSize',15);
                elseif (i >=10 && i<13)
                    ylabel('a2: Linear acceleration [m/s^2]','FontSize',15);
                end

                axis tight;
                grid on;
            end
        end
        fprintf('\n');
    end
end

fprintf('---------\n');
fprintf('Check done!\n');

