clc; close all; clear;

%% testOptions
Section1 = false;            % section for creating of fake subject;    

Section2 = true;            % section STEP1;  
plotTauComparison = true;

Section3 = true;             % section STEP2;

%% %%%%%%%%%%%%%%%%%%%%%%%%% ANALYSIS ASSUMPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%%
% COnsidering only a subjectList = [2,3,13,14];
% This analysis can be done only if:
%
% - SUBJECT 2       : URDF of the subject and measurements y_2;
% - SUBJECT 3       : URDF of subject+5kg (5kg positioned on the COM of 
%                     torso with zero Ixx, Iyy, Izz) and measurements y_3;
% - SUBJECT 13(fake): URDF of subject and measurements y_3;
% - SUBJECT 14(fake): URDF of 5kg  considering 0 mass on foot and leg and 
%                     5kg mass on torso (5kg positioned on the COM of 
%                     torso with zero Ixx, Iyy, Izz) with measurements y_3;
% - SUBJECT 15(fake): URDF of 5kg  considering 0 mass on foot and leg and 
%                     5kg mass on torso (5kg positioned on the COM of 
%                     torso with zero Ixx, Iyy, Izz) with measurements y_2;

%% %%%%%%%%%%%%%%%%%%%% CREATION OF FAKE URDF+MEAS %%%%%%%%%%%%%%%%%%%%%%%%
%%
if(Section1)
    creatingFakeSubject
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% GT ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add folders

dataFolder = './experiments/humanFixedBase/data/';
if (exist (dataFolder,'dir')==7)
    addpath(genpath(dataFolder));
end

helperFunctionsFolder = './experiments/humanFixedBase/helperFunctions/';
if (exist (helperFunctionsFolder,'dir')==7)
    addpath(genpath(helperFunctionsFolder));
end

intermediateDataFilesFolder = './experiments/humanFixedBase/intermediateDataFiles/';
if (exist(intermediateDataFilesFolder,'dir')==7)
    addpath(genpath(intermediateDataFilesFolder));
end

testsFolder = './experiments/humanFixedBase/tests/';
if (exist(testsFolder,'dir')==7)
    addpath(genpath(testsFolder));
end

%% folder for plots
figFolder = './experiments/humanFixedBase/plots';
if(exist(figFolder,'dir')==0)
    mkdir(figFolder);
end

%% add files .mat

fprintf('\nStarting synchronisedData computation\n');
synchroniseCaptureData
offset = 3.7852 *(pi/180);
synchronisedData(13,2).q(:,2) = synchronisedData(13,2).q(:,2) - offset*ones(size(synchronisedData(13,2).q(:,2)));
synchronisedData(3,2).q(:,2) = synchronisedData(3,2).q(:,2) - offset*ones(size(synchronisedData(3,2).q(:,2)));
synchronisedData(14,2).q(:,2) = synchronisedData(14,2).q(:,2) - offset*ones(size(synchronisedData(14,2).q(:,2)));
save('./experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat','synchronisedData');

% fprintf('\nStarting computeSubjectSpecificURDFPrams computation\n');
% computeSubjectSpecificURDFParams
% save('./experiments/humanFixedBase/data/subjectSizeParams.mat','subjectParams');
% % % % % 
% % % % % fprintf('\nStarting createUrdfModelFromSubjectParam computation\n'); 
% % % % % createUrdfModelFromSubjectParams

fprintf('\nStarting loadModelFromURDF computation\n');  
loadModelFromURDF

fprintf('\nStarting computeLinkSensorFrames computation\n');
computeLinkSensorFrames

fprintf('\nStarting organiseBERDYCompatibleSensorData computation\n');   
organiseBERDYCompatibleSensorData

%% MAP computation

fprintf('\nStarting MAP computation\n');
% selected subjects and trials
subjectList = [2,3,13,14,15];
trialList = 2;  

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d ',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d\n',trialID);

        %% load BERDY compatible data and sensor link transforms
        load('./experiments/humanFixedBase/intermediateDataFiles/BERDYFormattedSensorData.mat');

        currentTrial = BERDYFormattedSensorData(subjectID,trialID); 
        data = currentTrial.data;
        dataTime = currentTrial.data.dataTime;
        q = currentTrial.data.q';
        dq = currentTrial.data.dq';
        ddq = currentTrial.data.ddq';
        len = length(dataTime);

        load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');

        currentTrialSens = sensorLinkTransforms(subjectID,trialID);
        X_imu_2 = currentTrialSens.X_imu_2;
        XStar_fp_0 = currentTrialSens.XStar_fp_0;
        XStar_0_1 = currentTrialSens.XStar_0_1;

        %% build model valid for MAP case 
        load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');

        currentModel = humanThreeLinkModelFromURDF(subjectID);

        humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
        humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 

        modelSigma = 1e+5;
        dmodel  = currentModel.dmodel;                      % deterministic model
        dmodel  = autoTreeStochastic(dmodel, modelSigma);         % probabilistic model for D equation (added Sv and Sw)
        myModel = model(dmodel);

        %% Build sensor model    
        sens.parts    = {'leg'         ,'torso'};       %force of the forceplate is trasmitted into the leg
        sens.labels   = {'fts'         ,'imu'  };
        sens.ndof     = {6             ,6      };
        label_to_plot = {'fts'         ,'imu'  };
        ymodel_4sens  = humanThreeLinkSens(dmodel, sens);  
        ymodel_4sens  = humanThreeLinkSensStochastic(ymodel_4sens);     % probabilistic model for Y(q,dq) d = y (added Sy)
        mySens_4sens  = sensors(ymodel_4sens);  
        myMAP_4sens   = MAP(myModel, mySens_4sens);

        %% Build data.y anda data.Sy 
        % data.y are ordered in:  
        % - angular-linear notation 
        % - the form [f1 a2 ftx1 ftx2 ddq1 ddq2]
        % - sensor frame

        %===== data.y
        data.y_4sens  = [];
        for i = 1 : length(sens.labels)
             eval(['data.y_4sens  = [data.y_4sens  data.ys_sensFrame_' sens.labels{i} '];']);
        end
        % Add null external forces ftx = 0
        data.y_4sens  = [data.y_4sens, zeros(len,6*dmodel.NB)];
        % Add ddq measurements
        data.y_4sens  = [data.y_4sens, 0.*data.ddq];
        data.y_4sens = data.y_4sens';

        %===== data.Sy
        data.Sy_4sens = [];
        for i = 1 : length(myMAP_4sens.IDsens.sensorsParams.labels)
             data.Sy_4sens = [data.Sy_4sens; diag(myMAP_4sens.IDsens.sensorsParams.Sy{i})];
        end
        data.Sy_4sens = repmat(data.Sy_4sens, 1, len-1);
        data.Sy_4sens = [data.Sy_4sens data.Sy_4sens(:,end)];

        %% Build Ymatrix manually
        % Ymatrix has to be consistent with measurements form [f1 a2 ftx1 ftx2 ddq1 ddq2]

        Y_4sens = cell2mat(ymodel_4sens.Y);
        %Y_4sens = zeros (ymodel.m,26*dmodel.NB);
        Y_4sens(7:12,27:32) = X_imu_2;
        Y_4sens(13:18,20:25) = eye(6);
        Y_4sens(19:24,46:51) = eye(6);
        Y_4sens(25,26) = eye(1);
        Y_4sens(26,52) = eye(1);

        Ymatrix_4sens = cell(len,1);
        for i = 1 : len
            %the only row in Ymatrix that is time varying
            Y_4sens(1:6,13:18) = XStar_fp_0 * XStar_0_1{i};
            Ymatrix_4sens{i} = Y_4sens; 
        end

         %% for computing v
        ymodel_RNEA  = autoSensRNEA(dmodel);
        mySens_RNEA  = sensors(ymodel_RNEA);
        myRNEA       = RNEA(myModel, mySens_RNEA);

        y_RNEA_f = zeros(6*dmodel.NB, len);
        y_RNEA_ddq = zeros(dmodel.NB, len);
        fx = cell(dmodel.NB);

        %Ordering y_RNEA in the form [fx1 fx2 ddq1 ddq2]
        for i = 1 : dmodel.NB
            for j = 1 : len
                fx{i,1} = zeros(6,1); 
                y_RNEA_f(((1:6)+(i-1)*6), j) = [fx{i,1}];
                y_RNEA_ddq(i, j) = [ddq(i,j)];
            end
            y_RNEA = [y_RNEA_f ; 0.*y_RNEA_ddq];
        end

        d = zeros (26*myRNEA.IDmodel.modelParams.NB,len);
        v_RNEA = cell(len,1);
        resRNEA.tau = zeros(length(dataTime),dmodel.NB);

        %tau_cla = zeros(length(dataTime),dmodel.NB);
        for i = 1 : len
             myRNEA = myRNEA.setState(q(:,i),0.*dq(:,i));
             myRNEA = myRNEA.setY(y_RNEA(:,i));
             myRNEA = myRNEA.solveID();

             d(:,i) = myRNEA.d;

             v_RNEA{i,:} = myRNEA.v; 
             resRNEA.tau(i,:) = myRNEA.tau;
             
%              deltaL1 = 0.79699 * sin(q(1,i));
%              deltaL2 = 0.26923 * sin(q(1,i) + q(2,i));
%              tau_cla(i,:) = [ ((deltaL2+deltaL1) * 5 * -9.81), (deltaL2 * 5 * -9.81)];
        end 
        clear fx;

        %% Build bias b_Y manually
        % b_Y has to be consistent with Ymatrix

        load('./experiments/humanFixedBase/data/subjectSizeParams_fake.mat');

        currentParams = subjectParams(subjectID);

        footMass =  currentParams.footMass;
        posP_0 = [0; 0; (0.5*currentParams.footHeight)];
        footIxx =  currentParams.footIxx;
        footIyy =  currentParams.footIyy;
        footIzz =  currentParams.footIzz;

        b_Y_4sens = zeros (size(data.y_4sens)); 
        R_imu_2 = X_imu_2(1:3,1:3);

        a_G_grav = [0;0;0;0;0;-9.8100]; %Featherstone-like notation, in global reference
        X_G_0 = currentTrialSens.X_G_0;
        X_0_G = InverseAdjTransform(X_G_0);
        a_0_grav = X_0_G * a_G_grav;

        I_0 = createSpatialInertia(footIxx,footIyy,footIzz,footMass,posP_0);

        b_Y_4sens(1:6,1:len)   = repmat((-XStar_fp_0 * I_0 * a_0_grav),1,len);

        v = cell(size(q))';
        fx = zeros (6,1);
        fext    = cell(1,2);
        for i = 1 : dmodel.NB
             fext{i}    = fx;
        end

        % exploiting velocitiy v_RNEA coming from RNEA class computation
        for i = 1 : len 
            A =R_imu_2*v_RNEA{i,1}(1:3,2);
            B =((X_imu_2(4:6,1:3)*v_RNEA{i,1}(1:3,2))+(R_imu_2*v_RNEA{i,1}(4:6,2)));
            b_Y_4sens(10:12,i) = cross(A,B);
        end 

        clear A;
        clear B;
        clear footIxx;
        clear footIyy;
        clear footIzz;

        %% Computing MAP method
        fprintf('\nMAP computation with all sensors\n');    
        resMAP_4sens.d  = zeros(26*dmodel.NB,len);
        resMAP_4sens.Sd = cell(len,1);
        resMAP_4sens.Ymatrix = cell(len,1);
        resMAP_4sens.b_Y = zeros (myMAP_4sens.IDsens.sensorsParams.m,len);

        for i = 1 : len

            myMAP_4sens = myMAP_4sens.setState(data.q(i,:)', 0.*data.dq(i,:)');
            myMAP_4sens = myMAP_4sens.setY(data.y_4sens(:,i));
            myMAP_4sens = myMAP_4sens.setYmatrix(Ymatrix_4sens{i});
            myMAP_4sens = myMAP_4sens.setBias(b_Y_4sens(:,i));
            myMAP_4sens = myMAP_4sens.solveID();

            resMAP_4sens.d(:,i)       = myMAP_4sens.d;
            resMAP_4sens.Sd{i,1}      = myMAP_4sens.Sd;
            resMAP_4sens.Ymatrix{i,1} = myMAP_4sens.IDsens.sensorsParams.Y; 
            resMAP_4sens.b_Y(:,i)     = myMAP_4sens.IDsens.sensorsParams.bias;

             if mod(i-1,100) == 0
                    fprintf('Processing %d %% of the dataset\n', round(i/len*100));
             end
        end
        % ========end MAP
        
        %% Rearrange solution

        for i = 1 : dmodel.NB

            link = strrep(myMAP_4sens.IDmodel.modelParams.linkname{i}, '+', '_');
            joint = strrep(myMAP_4sens.IDmodel.modelParams.jointname{i}, '+', '_');

            % initialize variables
            eval(['resMAP_4sens.a_'    link ' = zeros(6, len );']);
            eval(['resMAP_4sens.Sa_'    link ' = cell(len, 1);']);

            eval(['resMAP_4sens.fB_'    link ' = zeros(6, len );']);
            eval(['resMAP_4sens.SfB_'    link ' = cell(len, 1);']);

            eval(['resMAP_4sens.f_'    link ' = zeros(6, len );']);
            eval(['resMAP_4sens.Sf_'     link ' = cell(len, 1);']);

            eval(['resMAP_4sens.tau_'    joint ' = zeros(1, len );']);
            eval(['resMAP_4sens.Stau_'    joint ' = zeros(len, 1);']);

            eval(['resMAP_4sens.fx_'    link ' = zeros(6, len );']);
            eval(['resMAP_4sens.Sfx_'    link ' = cell(len, 1);']);

            eval(['resMAP_4sens.d2q_'    joint ' = zeros(1, len );']);
            eval(['resMAP_4sens.Sd2q_'    joint ' = zeros(len, 1);']);

            for j = 1 : len
                 %a
                 ind  = '1 + 26*(i-1) : 26*(i-1) +  6';
                 eval(['resMAP_4sens.a_'    link '(:,j) = resMAP_4sens.d      (' ind '       ,j);'])
                 eval(['resMAP_4sens.Sa_'   link '{j,1} = resMAP_4sens.Sd{j,1}(' ind ',' ind ' );'])
                 %fB
                 ind  = '7 + 26*(i-1) : 26*(i-1) + 12';
                 eval(['resMAP_4sens.fB_'   link '(:,j) = resMAP_4sens.d      (' ind '        ,j);'])
                 eval(['resMAP_4sens.SfB_'  link '{j,1} = resMAP_4sens.Sd{j,1}(' ind ',' ind '  );'])
                 %f
                 ind  = '13 + 26*(i-1) : 26*(i-1) + 18';
                 eval(['resMAP_4sens.f_'    link '(:,j) = resMAP_4sens.d      (' ind '        ,j);'])
                 eval(['resMAP_4sens.Sf_'   link '{j,1} = resMAP_4sens.Sd{j,1}(' ind ',' ind '  );'])
                 %tau
                 ind  = '19 + 26*(i-1) : 26*(i-1) + 19';
                 eval(['resMAP_4sens.tau_'  joint '(:,j) = resMAP_4sens.d      (' ind '        ,j);'])
                 eval(['resMAP_4sens.Stau_' joint '(j,1) = resMAP_4sens.Sd{j,1}(' ind ',' ind '  );'])
                 %fx
                 ind  = '20 + 26*(i-1) : 26*(i-1) + 25';
                 eval(['resMAP_4sens.fx_'   link '(:,j) = resMAP_4sens.d      (' ind '        ,j);'])
                 eval(['resMAP_4sens.Sfx_'  link '{j,1} = resMAP_4sens.Sd{j,1}(' ind ',' ind '  );'])
                 %ddq
                 ind  = '26 + 26*(i-1) : 26*(i-1) + 26';
                 eval(['resMAP_4sens.d2q_'  joint '(:,j) = resMAP_4sens.d      (' ind '        ,j);'])
                 eval(['resMAP_4sens.Sd2q_' joint '(j,1) = resMAP_4sens.Sd{j,1}(' ind ',' ind '  );'])   
             end
        end 

         %% Organising into a structure    

        GTresults(subjectID,trialID).dataTime = dataTime;
        GTresults(subjectID,trialID).resMAP_4sens = resMAP_4sens;
        GTresults(subjectID,trialID).resRNEA.tau = resRNEA.tau;
        GTresults(subjectID,trialID).resMAP_4sens.data.y = data.y_4sens;
        clear res;

     end
fprintf('\n');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validation of the following equation:
% tau(subject+5kg,y_subject+5kg) - tau(subject,y_subject+5kg) = tauHat(5kg)
%
% tauHat(5kg) obtained from this equation has to be identical to tau(5kg)
% coming from subject14.  
%
% Note1: this relation will be verified only if the center of mass doesn't
% change in URDFs.
% Note2: this relation has to be identical with and without q,dq,ddq=0;
%%
if(Section2)
    %% joint1
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    title('\tau_{5kg} RNEA for joint 1');

    dataTimeMin = GTresults(13,2).dataTime;

    % COMPUTING TAU HAT 5KG --> RNEA
    tau1_ankleSubj3_RNEA  = GTresults(3,2).resRNEA.tau(:,1);
    tau1_ankleSubj13_RNEA = GTresults(13,2).resRNEA.tau(:,1);

    tau1Hat_5kg_RNEA = tau1_ankleSubj3_RNEA' - tau1_ankleSubj13_RNEA';

    % COMPUTING TAU 5KG --> RNEA
    tau1_5kg_RNEA = GTresults(14,2).resRNEA.tau(:,1);

    plot1 = plot(dataTimeMin,tau1Hat_5kg_RNEA,'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2= plot(dataTimeMin,tau1_5kg_RNEA,'--','lineWidth',1.5); hold on;
    set(plot2,'color',[0 0 1]);

    leg = legend('$\hat\tau_{5kg}$','$\tau_{5kg}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',30);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_1$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;

    if(plotTauComparison)
       save2pdf(fullfile(figFolder,'joint1_tauRNEAcomparison'),fig,600);
    end

    %% joint2

    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    title('\tau_{5kg} RNEA for joint 2');

    dataTimeMin = GTresults(13,2).dataTime;

    % COMPUTING TAU HAT 5KG --> RNEA
    tau2_hipSubj3_RNEA  = GTresults(3,2).resRNEA.tau(:,2);
    tau2_hipSubj13_RNEA = GTresults(13,2).resRNEA.tau(:,2);

    tau2Hat_5kg_RNEA = tau2_hipSubj3_RNEA - tau2_hipSubj13_RNEA;

    % COMPUTING TAU 5KG
    tau2_5kg_RNEA = GTresults(14,2).resRNEA.tau(:,2);
    tau2_5kg_RNEA_y2 = GTresults(15,2).resRNEA.tau(:,2);

    plot1 = plot(dataTimeMin,tau2Hat_5kg_RNEA,'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2= plot(dataTimeMin,tau2_5kg_RNEA,'--','lineWidth',1.5); hold on;
    set(plot2,'color',[0 0 1]);

    leg = legend('$\hat\tau_{5kg}$','$\tau_{5kg}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',30);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
     
    if(plotTauComparison)
       save2pdf(fullfile(figFolder, 'joint2_tauRNEAcomparison'),fig,600);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if(Section3)
%% joint1
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
 
    dataT = synchronisedData(13,2).dataTime;
    
    tau1_ankleSubj2_RNEA = GTresults(2,2).resRNEA.tau(:,1);
    
    subplot(311)
    plot1 = plot(dataT,tau1_ankleSubj13_RNEA(1:length(dataT)),'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2= plot(dataT,tau1_ankleSubj2_RNEA(1:length(dataT)),'lineWidth',1.5); hold on;
    set(plot2,'color',[0 0 1]);
    plot3= plot(dataT,tau1_5kg_RNEA(1:length(dataT)),'lineWidth',1.5); hold on;
    set(plot3,'color',[0 1 0]);
    title(sprintf('Comparison with model variance %d',modelSigma));
    
    leg = legend('$\tau_{13}$','$\tau_{2}$','$\tau_{5kg}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_1$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
   
    
    q1Subj2 =(synchronisedData(2, 2).q(:,1));
    q1Subj2Cut = q1Subj2(1:length(dataT));
    
    q1Subj13 =synchronisedData(13, 2).q(:,1);
    
    
    subplot(312)
    tau1_ankleSubj2_RNEA = tau1_ankleSubj2_RNEA(1:length(dataT));
    diffTau1_RNEA = tau1_ankleSubj13_RNEA-tau1_ankleSubj2_RNEA-tau1_5kg_RNEA;
    
    tau1_ankleSubj13_MAP  = GTresults(13,2).resMAP_4sens.tau_ankle;
    tau1_ankleSubj2_MAP = GTresults(2, 2).resMAP_4sens.tau_ankle;
    tau1_ankleSubj2_MAP  = tau1_ankleSubj2_MAP(1:length(dataT));
    
    diffTau1_MAP = tau1_ankleSubj13_MAP-tau1_ankleSubj2_MAP-tau1_5kg_RNEA';
  
%     integral1_RNEA = trapz(abs(diffTau1_RNEA))
%     integral1_MAP  = trapz(abs(diffTau1_MAP))

    
    plot1 = plot(dataT,diffTau1_RNEA(1:length(dataT)),'k','lineWidth',1.5); hold on;
%     set(plot1,'color',[1 0 0]);
    plot2= plot(dataT,diffTau1_MAP(1:length(dataT)),'m','lineWidth',1.5); hold on;
%     set(plot2,'color',[0 0 1]);
    
    leg = legend('$Diff_{RNEA}$','$Diff_{MAP}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_1$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    

    subplot(313)
    plot1 = plot(dataT,(180/pi)*q1Subj2Cut,'lineWidth',1.5); hold on;
    set(plot1,'color',[0 0 1]);
    plot2= plot(dataT,(180/pi)*q1Subj13,'lineWidth',1.5); hold on;
    set(plot2,'color',[1 0 0]);
    
    leg = legend('$q_{2}$','$q_{13}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$q_1$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
%% joint2
    
fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
 
    dataT = synchronisedData(13,2).dataTime;
    
    tau2_hipSubj2_RNEA = GTresults(2,2).resRNEA.tau(:,2);
    
    subplot(311)
    plot1 = plot(dataT,tau2_hipSubj13_RNEA(1:length(dataT)),'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2= plot(dataT,tau2_hipSubj2_RNEA(1:length(dataT)),'lineWidth',1.5); hold on;
    set(plot2,'color',[0 0 1]);
    plot3= plot(dataT,tau2_5kg_RNEA(1:length(dataT)),'lineWidth',1.5); hold on;
    set(plot3,'color',[0 1 0]);
    title(sprintf('Comparison with model variance %d',modelSigma));
    
    leg = legend('$\tau_{13}$','$\tau_{2}$','$\tau_{5kg}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
   
    
    xxx=(synchronisedData(2, 2).q(:,2));
    xxx = xxx(1:length(dataT));
    
    subplot(312)
    tau2_hipSubj2_RNEA = tau2_hipSubj2_RNEA(1:length(dataT));
    diffTau2_RNEA = tau2_hipSubj13_RNEA-tau2_hipSubj2_RNEA-tau2_5kg_RNEA;
    
    tau2_hipSubj13_MAP  = GTresults(13,2).resMAP_4sens.tau_hip;
    tau2_hipSubj2_MAP = GTresults(2, 2).resMAP_4sens.tau_hip;
    tau2_hipSubj2_MAP  = tau2_hipSubj2_MAP(1:length(dataT));
    
    diffTau2_MAP = tau2_hipSubj13_MAP-tau2_hipSubj2_MAP-tau2_5kg_RNEA';
  
    integral2_RNEA = trapz(abs(diffTau2_RNEA))
    integral2_MAP  = trapz(abs(diffTau2_MAP))
    
    plot1 = plot(dataT,diffTau2_RNEA(1:length(dataT)),'k','lineWidth',1.5); hold on;
    %set(plot1,'color',[1 0 0]);
    plot2= plot(dataT,diffTau2_MAP(1:length(dataT)),'m','lineWidth',1.5); hold on;
    %set(plot2,'color',[0 0 1]);
    
    leg = legend('$Diff_{RNEA}$','$Diff_{MAP}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    
    subplot(313)
    plot1 = plot(dataT,(180/pi)*xxx,'lineWidth',1.5); hold on;
    set(plot1,'color',[0 0 1]);
    plot2= plot(dataT,(180/pi)*synchronisedData(13, 2).q(:,2),'lineWidth',1.5); hold on;
    set(plot2,'color',[1 0 0]);
    
    leg = legend('$q_{2}$','$q_{13}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    %% plotting Tau2_RNEA in (q1+q2)
    
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    %title('\tau_2 in (q_1+q_2)');
    
    % q1+q2 for Subject2
    q1Subj2 =(synchronisedData(2, 2).q(:,1));
    q1Subj2Cut = q1Subj2(1:length(dataT));
    q2Subj2 =(synchronisedData(2, 2).q(:,2));
    q2Subj2Cut = q2Subj2(1:length(dataT));
    
    qSumSubj2 = q1Subj2Cut + q2Subj2Cut;
    qSumSubj2 = (180/pi)* qSumSubj2;
    
    % q1+q2 for Subject13
    q1Subj13 =(synchronisedData(13, 2).q(:,1));
    q1Subj13Cut = q1Subj13(1:length(dataT));
    q2Subj13 =(synchronisedData(13, 2).q(:,2));
    q2Subj13Cut = q2Subj13(1:length(dataT));
    
    qSumSubj13 = q1Subj13Cut + q2Subj13Cut;
    qSumSubj13 = (180/pi)* qSumSubj13;
    
    
    
    subplot(211)
    plot1 = plot(dataT,tau2_hipSubj13_RNEA(1:length(dataT)),'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2= plot(dataT,tau2_hipSubj2_RNEA(1:length(dataT)),'lineWidth',1.5); hold on;
    set(plot2,'color',[0 0 1]);
    title('RNEA computation');
    
    leg = legend('$\tau_{13}$','$\tau_{2}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;

    
    subplot(212)
    
    coeff_tau2_hipSubj2_RNEA = polyfit(qSumSubj2,tau2_hipSubj2_RNEA(1:length(dataT)), 4);
    fit_tau2_hipSubj2_RNEA = polyval(coeff_tau2_hipSubj2_RNEA, qSumSubj2); 
    
    coeff_tau2_hipSubj13_RNEA = polyfit(qSumSubj13,tau2_hipSubj13_RNEA(1:length(dataT)), 4);
    fit_tau2_hipSubj13_RNEA = polyval(coeff_tau2_hipSubj13_RNEA, qSumSubj13);
    
    plot1 = plot(qSumSubj13,tau2_hipSubj13_RNEA(1:length(dataT)),'.','lineWidth',1.5); hold on;
    set(plot1,'color',[1 0.400000005960464 0.400000005960464]);
    plot2 = plot(qSumSubj13,fit_tau2_hipSubj13_RNEA,'lineWidth',3); hold on;
    set(plot2,'color',[1 0 0]);
    
    plot3 = plot(qSumSubj2,tau2_hipSubj2_RNEA(1:length(dataT)),'.','lineWidth',0.1); hold on;
    set(plot3,'color',[0.301960796117783 0.745098054409027 0.933333337306976]);
    plot4 = plot(qSumSubj2,fit_tau2_hipSubj2_RNEA,'lineWidth',3); hold on;
    set(plot4,'color',[0 0 1]);
    

    leg = legend('$\tau_{13}$','$\tau_{fit13}$','$\tau_{2}$','$\tau_{fit2}$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('$q_1 + q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    % computing the difference
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    
    qs = 2:0.01:32;
    fit2RNEA_qs = polyval(coeff_tau2_hipSubj2_RNEA, qs);
    fit13RNEA_qs =  polyval(coeff_tau2_hipSubj13_RNEA, qs);
    
    subplot(211)
    plot1 = plot(qs,fit13RNEA_qs,'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2= plot(qs,fit2RNEA_qs,'lineWidth',1.5); hold on;
    set(plot2,'color',[0 0 1]);
    title('RNEA computation')
    
    leg = legend('$\tau_{fit13}$','$\tau_{fit2}$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('$q_1 + q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    
    diffInterpRNEA = fit13RNEA_qs - fit2RNEA_qs;
    
    subplot(212)
    plot1 = plot(qs,diffInterpRNEA,'lineWidth',1.5); hold on;
    set(plot1,'color',[0 1 0]);

    leg = legend('$diff$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('$q_1 + q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;

    %% plotting Tau2_MAP in (q1+q2)
    
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    %title('\tau_2 in (q_1+q_2)');
    
    % q1+q2 for Subject2
    q1Subj2 =(synchronisedData(2, 2).q(:,1));
    q1Subj2Cut = q1Subj2(1:length(dataT));
    q2Subj2 =(synchronisedData(2, 2).q(:,2));
    q2Subj2Cut = q2Subj2(1:length(dataT));
    
    qSumSubj2 = q1Subj2Cut + q2Subj2Cut;
    qSumSubj2 = (180/pi)* qSumSubj2;
    
    % q1+q2 for Subject13
    q1Subj13 =(synchronisedData(13, 2).q(:,1));
    q1Subj13Cut = q1Subj13(1:length(dataT));
    q2Subj13 =(synchronisedData(13, 2).q(:,2));
    q2Subj13Cut = q2Subj13(1:length(dataT));
    
    qSumSubj13 = q1Subj13Cut + q2Subj13Cut;
    qSumSubj13 = (180/pi)* qSumSubj13;
    
    tau2_hipSubj13_MAP  = GTresults(13,2).resMAP_4sens.tau_hip;
    tau2_hipSubj2_MAP = GTresults(2, 2).resMAP_4sens.tau_hip;
    tau2_hipSubj2_MAP  = tau2_hipSubj2_MAP(1:length(dataT));
    
    
    subplot(211)
    plot1 = plot(dataT,tau2_hipSubj13_MAP,'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2= plot(dataT,tau2_hipSubj2_MAP,'lineWidth',1.5); hold on;
    set(plot2,'color',[0 0 1]);
    title(sprintf('model variance = %d, MAP computation',modelSigma))
    
    leg = legend('$\tau_{13}$','$\tau_{2}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;

    
    subplot(212)
    coeff_tau2_hipSubj2_MAP = polyfit(qSumSubj2',tau2_hipSubj2_MAP, 7);
    fit_tau2_hipSubj2_MAP = polyval(coeff_tau2_hipSubj2_MAP, qSumSubj2); 
    
    coeff_tau2_hipSubj13_MAP = polyfit(qSumSubj13',tau2_hipSubj13_MAP, 7);
    fit_tau2_hipSubj13_MAP = polyval(coeff_tau2_hipSubj13_MAP, qSumSubj13);
    
    plot1 = plot(qSumSubj13,tau2_hipSubj13_MAP,'.','lineWidth',1.5); hold on;
    set(plot1,'color',[1 0.400000005960464 0.400000005960464]);
    plot2 = plot(qSumSubj13,fit_tau2_hipSubj13_MAP,'lineWidth',3); hold on;
    set(plot2,'color',[1 0 0]);
    
    plot3 = plot(qSumSubj2,tau2_hipSubj2_MAP,'.','lineWidth',0.1); hold on;
    set(plot3,'color',[0.301960796117783 0.745098054409027 0.933333337306976]);
    plot4 = plot(qSumSubj2,fit_tau2_hipSubj2_MAP,'lineWidth',3); hold on;
    set(plot4,'color',[0 0 1]);
    
    leg = legend('$\tau_{13}$','$\tau_{fit13}$','$\tau_{2}$','$\tau_{fit2}$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('$q_1 + q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    % computing the difference
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    
    qs = 2:0.01:32;
    fit2MAP_qs = polyval(coeff_tau2_hipSubj2_MAP, qs);
    fit13MAP_qs =  polyval(coeff_tau2_hipSubj13_MAP, qs);
    
    subplot(211)
    plot1 = plot(qs,fit13MAP_qs,'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2= plot(qs,fit2MAP_qs,'lineWidth',1.5); hold on;
    set(plot2,'color',[0 0 1]);
    title('MAP computation')
    
    leg = legend('$\tau_{fit13}$','$\tau_{fit2}$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('$q_1 + q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    
    diffInterpMAP = fit13MAP_qs - fit2MAP_qs;
    
    subplot(212)
    plot1 = plot(qs,diffInterpMAP,'lineWidth',1.5); hold on;
    set(plot1,'color',[0 1 0]);

    
    leg = legend('$diff$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('$q_1 + q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    %% tau_5kg

    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    
    q1Subj14 =(synchronisedData(14, 2).q(:,1));
    q2Subj14 =(synchronisedData(14, 2).q(:,2));
    
    qSumSubj14 = q1Subj14 + q2Subj14;
    qSumSubj14 = (180/pi)* qSumSubj14;
    
    tau2_5kg_RNEA = GTresults(14,2).resRNEA.tau(:,2);
    
    subplot (311)
    plot1 = plot(dataTimeMin,tau2Hat_5kg_RNEA,'lineWidth',1.5); hold on;
    set(plot1,'color',[0.929411768913269 0.694117665290833 0.125490203499794]);
    leg = legend('$\tau_{5kg}$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    title('RNEA of 5kg');
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;

    
    subplot (312)
    plot2 = plot(qSumSubj14,tau2_5kg_RNEA);
    set(plot2,'color',[0.929411768913269 0.694117665290833 0.125490203499794]);
    leg = legend('$\tau_{5kg}$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('$q1+q2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;

    coeff_tau2_5kg_RNEA = polyfit(qSumSubj14,tau2_5kg_RNEA, 7);
    fit_tau2_5kg_RNEA = polyval(coeff_tau2_5kg_RNEA, qSumSubj14);
    qs = 2:0.01:32;
    fit14RNEA_qs = polyval(coeff_tau2_5kg_RNEA, qs);

    
    subplot (313)
    plot3 = plot(qs,fit14RNEA_qs);
    set(plot3,'color',[0.929411768913269 0.694117665290833 0.125490203499794]);
    leg = legend('$\tau_{5kg,fit}$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',20);
    xlabel('$q1+q2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    
    %%  SUM Angles
    
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    
    q1Sub13 = synchronisedData(13, 2).q(:,1)  * (180/pi);
    q2Sub13 = synchronisedData(13, 2).q(:,2) * (180/pi);
    
    q1Sub2 = synchronisedData(2, 2).q(:,1);
    q1Sub2 = q1Sub2(1:length(q1Sub13)) * (180/pi);
    q2Sub2 = synchronisedData(2, 2).q(:,2);
    q2Sub2 = q2Sub2(1:length(q1Sub13)) * (180/pi);
    

    subplot (311)
    plot1 = plot(dataTimeMin,q1Sub2,'lineWidth',1.5); hold on;
    set(plot1,'color',[0 0 1]);
    plot2 = plot(dataTimeMin,q1Sub13,'lineWidth',1.5); hold on;
    set(plot2,'color',[1 0 0]);
    leg = legend('subj2','subj13','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',15);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$q_1$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    
    subplot (312)
    plot1 = plot(dataTimeMin,q2Sub2,'lineWidth',1.5); hold on;
    set(plot1,'color',[0 0 1]);
    plot2 = plot(dataTimeMin,q2Sub13,'lineWidth',1.5); hold on;
    set(plot2,'color',[1 0 0]);
%     leg = legend('subj2','subj13','Location','northeast');
%     set(leg,'Interpreter','latex');
    set(leg,'FontSize',15);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    

    subplot (313)
    plot1 = plot(dataTimeMin,q1Sub13 + q2Sub13,'lineWidth',1.5); hold on;
    set(plot1,'color',[0 0 1]);
    plot2 = plot(dataTimeMin,q1Sub2 + q2Sub2,'lineWidth',1.5); hold on;
    set(plot2,'color',[1 0 0]);
    set(leg,'FontSize',15);
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$q_1 + q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    %% angles theta_IMU from rough data
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    %title('theta');
    
    % SUBJECT2
    theta_subj2 = -atan2(synchronisedData(2, 2).aLin_imu_imu(:,3), synchronisedData(2, 2).aLin_imu_imu(:,1));
    theta_subj2 = theta_subj2 * (180/pi);
    
    % SUBJECT13
    theta_subj13 = -atan2(synchronisedData(13, 2).aLin_imu_imu(:,3), synchronisedData(13, 2).aLin_imu_imu(:,1));
    theta_subj13 = theta_subj13 * (180/pi);
    
    subplot (311)
    plot1 = plot(dataTimeMin,theta_subj13,'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2 = plot(dataTimeMin,theta_subj2(1:length(theta_subj13)),'lineWidth',1.5);
    set(plot2,'color',[0 0 1]);
%     leg = legend('subj13','subj2','Location','southeast');
%     set(leg,'Interpreter','latex');
%     set(leg,'FontSize',15);
    title('\theta from IMU analysis');
     
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\theta_{IMU}$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    subplot (312)
    plot1 = plot(dataTimeMin,qSumSubj13,'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2 = plot(dataTimeMin,qSumSubj2(1:length(theta_subj13)),'lineWidth',1.5);
    set(plot2,'color',[0 0 1]);
    leg = legend('subj13','subj2','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',15);
    
    xlabel('Time [s]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$q_1 + q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    subplot (313)
    plot1 = plot(qSumSubj13,theta_subj13,'lineWidth',1.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2 = plot(qSumSubj2,theta_subj2(1:length(theta_subj13)),'lineWidth',1.5);
    set(plot2,'color',[0 0 1]);
%     leg = legend('subj13','subj2','Location','southeast');
%     set(leg,'Interpreter','latex');
%     set(leg,'FontSize',15);
    
    xlabel('$q_1 + q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\theta_{IMU}$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    axis equal;
    grid on;
    
    %% Final plot
    
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    %title('theta');

    plot1 = plot(qs,diffInterpMAP,'lineWidth',1.5); hold on;
    plot2 = plot(qs,diffInterpRNEA,'lineWidth',1.5); hold on;
    plot3 = plot(qs,fit14RNEA_qs,'lineWidth',1.5); hold on;

    leg = legend('MAP','RNEA','RNEA5kg','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',15);    
    
    xlabel('$q_1 + q_2$ [deg]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    ylabel('$\tau_{2}$ [Nm]','HorizontalAlignment','center',...
           'FontWeight','bold',...
           'FontSize',20,...
           'Interpreter','latex');
    axis tight;
    grid on;
    
    
end