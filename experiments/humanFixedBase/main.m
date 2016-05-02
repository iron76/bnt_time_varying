
clear;clc;close all;

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

%% add files .mat

fprintf('Starting synchronisedData computation\n');
file = './experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat';
if(exist(file,'file')==2)
    load(file);
else    
    synchroniseCaptureData
end
%
fprintf('\nStarting computeSubjectSpecificURDFPrams computation\n');
file = './experiments/humanFixedBase/data/subjectSizeParams.mat';
if(exist(file,'file')==2)
    load(file);
else    
    computeSubjectSpecificURDFParams
end
%
fprintf('\nStarting createUrdfModelFromSubjectParam computation\n');
folder = './human_models/';
if(exist(folder,'dir')==7)
    addpath(genpath(folder))
else    
    createUrdfModelFromSubjectParams
end
%
fprintf('\nStarting loadModelFromURDF computation\n');  
    loadModelFromURDF
%
fprintf('\nStarting computeLinkSensorFrames computation\n');
file = './experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat';
if(exist(file,'file')==2)
    load(file);
else    
    computeLinkSensorFrames
end
%
fprintf('\nStarting organiseBERDYCompatibleSensorData computation\n');
file = './experiments/humanFixedBase/intermediateDataFiles/BERDYFormattedSensorData.mat';
if(exist(file,'file')==2)
    load(file);
else    
    organiseBERDYCompatibleSensorData
end

%% selected subjects and trials
subjectList = 2:3;
trialList = 1;  

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d ',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d\n',trialID);
        
        %% load BERDY compatible data and sensor link transforms
        load('./experiments/humanFixedBase/intermediateDataFiles/BERDYFormattedSensorData.mat');

        currentTrial = BERDYFormattedSensorData(subjectID,trialID); 
        data = currentTrial.data;
        dataTime = data.dataTime;
        q = currentTrial.data.q';
        dq = currentTrial.data.dq';
        ddq = currentTrial.data.ddq';
        len = length(dataTime);

        load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');

        currentTrialSens = sensorLinkTransforms(subjectID,trialID);
        X_imu_2 = currentTrialSens.X_imu_2;
        XStar_fp_0 = currentTrialSens.XStar_fp_0;
        XStar_0_1 = currentTrialSens.XStar_0_1;
% 
%         % =====structure from files
%         data.parts    = {'leg'         ,'torso'};
%         data.labels   = {'fts'         ,'imu'  };
%         data.ndof     = {6             ,6      };
%         data.index    = {'1:6'         ,'1:6'  };

        %% build model valid for all MAP case 
        load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');

        currentModel = humanThreeLinkModelFromURDF(subjectID);

        humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
        humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 
        
        dmodel  = currentModel.dmodel;                      % deterministic model
        dmodel  = autoTreeStochastic(dmodel, 1e-4);         % probabilistic model for D equation (added Sv and Sw)
        myModel = model(dmodel);
        
        
        %% ======================= MAP with 2 sensors(ddq+ftx(1&2)) ===========================
        %% Build sensor model                  
        sens_2sens.parts    = {};      
        sens_2sens.labels   = {};
        sens_2sens.ndof     = {};
        ymodel_2sens  = humanThreeLinkSens(dmodel, sens_2sens);  
        ymodel_2sens  = humanThreeLinkSensStochastic(ymodel_2sens);     
        mySens_2sens  = sensors(ymodel_2sens);  
        myMAP_2sens   = MAP(myModel, mySens_2sens);
        
        %% Build data.y anda data.Sy 
        % data.y are ordered in:  
        % - angular-linear notation 
        % - the form [ftx1 ftx2 ddq1 ddq2]
        % - sensor frame

        %===== data.y
        data.y_2sens  = [];

        % Add null external forces ftx = 0
        data.y_2sens  = [data.y_2sens, zeros(len,6*dmodel.NB)];
        % Add ddq measurements
        data.y_2sens  = [data.y_2sens, data.ddq];
        data.y_2sens = data.y_2sens';

        %===== data.Sy
%         data.Sy = [];
%         data.Sy = repmat(data.Sy, 1, len-1);
%         data.Sy = [data.Sy data.Sy(:,end)];
       
        %% Build Ymatrix manually
        % Y has to be consistent with measurements form [ftx1 ftx2 ddq1 ddq2]

        Y_2sens = zeros(myMAP_2sens.IDsens.m, 26*dmodel.NB);
        Y_2sens(1:6,20:25) = eye(6);
        Y_2sens(7:12,46:51) = eye(6);
        Y_2sens(13,26) = eye(1);
        Y_2sens(14,52) = eye(1);
        
        %% Build bias manually
        % b_Y has to be consistent with Y

        b_Y_2sens = zeros(myMAP_2sens.IDsens.m,len);
        
        %% Computing MAP method
         fprintf('MAP computation with 2 sensors\n');
         resMAP_2sens.d  = zeros(26*dmodel.NB,len);
         resMAP_2sens.Sd = cell(len,1);
         resMAP_2sens.Ymatrix = cell(len,1);
         resMAP_2sens.b_Y = zeros (myMAP_2sens.IDsens.sensorsParams.m,len);
        
        for i = 1 : len

            myMAP_2sens = myMAP_2sens.setState(data.q(i,:)', data.dq(i,:)');
            myMAP_2sens = myMAP_2sens.setY(data.y_2sens(:,i));
            myMAP_2sens = myMAP_2sens.setYmatrix(Y_2sens);
            myMAP_2sens = myMAP_2sens.setBias(b_Y_2sens(:,i));
            myMAP_2sens = myMAP_2sens.solveID();

            resMAP_2sens.d(:,i)       = myMAP_2sens.d;
            resMAP_2sens.Sd{i,1}      = myMAP_2sens.Sd;
            resMAP_2sens.Ymatrix{i,1} = myMAP_2sens.IDsens.sensorsParams.Y; 
            resMAP_2sens.b_Y(:,i)     = myMAP_2sens.IDsens.sensorsParams.bias;
            
             if mod(i-1,100) == 0
                    fprintf('Processing %d %% of the dataset\n', round(i/len*100));
             end
        end
        % ========end MAP
        
       %% Rearrange solution

        for i = 1 : dmodel.NB
             
            link = strrep(myMAP_2sens.IDmodel.modelParams.linkname{i}, '+', '_');
            joint = strrep(myMAP_2sens.IDmodel.modelParams.jointname{i}, '+', '_');
             
            % initialize variables
            eval(['resMAP_2sens.a_'    link ' = zeros(6, len );']);
            eval(['resMAP_2sens.Sa_'    link ' = cell(len, 1);']);
            
            eval(['resMAP_2sens.fB_'    link ' = zeros(6, len );']);
            eval(['resMAP_2sens.SfB_'    link ' = cell(len, 1);']);
            
            eval(['resMAP_2sens.f_'    link ' = zeros(6, len );']);
            eval(['resMAP_2sens.Sf_'     link ' = cell(len, 1);']);
            
            eval(['resMAP_2sens.tau_'    joint ' = zeros(1, len );']);
            eval(['resMAP_2sens.Stau_'    joint ' = zeros(len, 1);']);
            
            eval(['resMAP_2sens.fx_'    link ' = zeros(6, len );']);
            eval(['resMAP_2sens.Sfx_'    link ' = cell(len, 1);']);
            
            eval(['resMAP_2sens.d2q_'    joint ' = zeros(1, len );']);
            eval(['resMAP_2sens.Sd2q_'    joint ' = zeros(len, 1);']);
             
            
            for j = 1 : len
                 %a
                 ind  = '1 + 26*(i-1) : 26*(i-1) +  6';
                 eval(['resMAP_2sens.a_'    link '(:,j) = resMAP_2sens.d      (' ind '       ,j);'])
                 eval(['resMAP_2sens.Sa_'   link '{j,1} = resMAP_2sens.Sd{j,1}(' ind ',' ind ' );'])
                 %fB
                 ind  = '7 + 26*(i-1) : 26*(i-1) + 12';
                 eval(['resMAP_2sens.fB_'   link '(:,j) = resMAP_2sens.d      (' ind '        ,j);'])
                 eval(['resMAP_2sens.SfB_'  link '{j,1} = resMAP_2sens.Sd{j,1}(' ind ',' ind '  );'])
                 %f
                 ind  = '13 + 26*(i-1) : 26*(i-1) + 18';
                 eval(['resMAP_2sens.f_'    link '(:,j) = resMAP_2sens.d      (' ind '        ,j);'])
                 eval(['resMAP_2sens.Sf_'   link '{j,1} = resMAP_2sens.Sd{j,1}(' ind ',' ind '  );'])
                 %tau
                 ind  = '19 + 26*(i-1) : 26*(i-1) + 19';
                 eval(['resMAP_2sens.tau_'  joint '(:,j) = resMAP_2sens.d      (' ind '        ,j);'])
                 eval(['resMAP_2sens.Stau_' joint '(j,1) = resMAP_2sens.Sd{j,1}(' ind ',' ind '  );'])
                 %fx
                 ind  = '20 + 26*(i-1) : 26*(i-1) + 25';
                 eval(['resMAP_2sens.fx_'   link '(:,j) = resMAP_2sens.d      (' ind '        ,j);'])
                 eval(['resMAP_2sens.Sfx_'  link '{j,1} = resMAP_2sens.Sd{j,1}(' ind ',' ind '  );'])
                 %ddq
                 ind  = '26 + 26*(i-1) : 26*(i-1) + 26';
                 eval(['resMAP_2sens.d2q_'  joint '(:,j) = resMAP_2sens.d      (' ind '        ,j);'])
                 eval(['resMAP_2sens.Sd2q_' joint '(j,1) = resMAP_2sens.Sd{j,1}(' ind ',' ind '  );'])   
             end
        end 
        

        %% ============================= MAP with 3 sensors(ddq+ftx(1&2)+fp) ================================
        %%
        %% Build sensor model    
        sens.parts    = {'leg'};       %force of the forceplate is trasmitted into the leg
        sens.labels   = {'fts'};
        sens.ndof     = {6};
        %label_to_plot = {'imu'};
        ymodel_3sens  = humanThreeLinkSens(dmodel, sens);  
        ymodel_3sens  = humanThreeLinkSensStochastic(ymodel_3sens);     % probabilistic model for Y(q,dq) d = y (added Sy)
        mySens_3sens  = sensors(ymodel_3sens);  
        myMAP_3sens   = MAP(myModel, mySens_3sens);

        %% Build data.y anda data.Sy 
        % data.y are ordered in:  
        % - angular-linear notation 
        % - the form [a2 ftx1 ftx2 ddq1 ddq2]
        % - sensor frame

        %===== data.y
        data.y_3sens  = [];
        for i = 1 : length(sens.labels)
             eval(['data.y_3sens  = [data.y_3sens  data.ys_sensFrame_' sens.labels{i} '];']);
        end
        % Add null external forces ftx = 0
        data.y_3sens  = [data.y_3sens, zeros(len,6*dmodel.NB)];
        % Add ddq measurements
        data.y_3sens  = [data.y_3sens, data.ddq];
        data.y_3sens = data.y_3sens';

% %         %===== data.Sy
% %         data.Sy_3sens = [];
% %         for i = 1 : length(myMAP_3sens.IDsens.sensorsParams.labels)
% %              data.Sy_3sens = [data.Sy_3sens; diag(myMAP_3sens.IDsens.sensorsParams.Sy{i})];
% %         end
% %         data.Sy_3sens = repmat(data.Sy_3sens, 1, len-1);
% %         data.Sy_3sens = [data.Sy_3sens data.Sy_3sens(:,end)];

        %% Build Ymatrix manually
        % Ymatrix has to be consistent with measurements form [a2 ftx1 ftx2 ddq1 ddq2]

        Y_3sens = cell2mat(ymodel_3sens.Y);
        %Y_3sens = zeros (ymodel.m,26*dmodel.NB);
        Y_3sens(7:12,20:25) = eye(6);
        Y_3sens(13:18,46:51) = eye(6);
        Y_3sens(19,26) = eye(1);
        Y_3sens(20,52) = eye(1);

        

        Ymatrix_3sens = cell(len,1);
        for i = 1 : len
            %the only row in Ymatrix that is time varying
            Y_3sens(1:6,13:18) = XStar_fp_0 * XStar_0_1{i};
            Ymatrix_3sens{i} = Y_3sens; 
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
            y_RNEA1 = [y_RNEA_f ; y_RNEA_ddq];
            y_RNEA2 = [y_RNEA_f ; y_RNEA_ddq];
        end

        d = zeros (26*myRNEA.IDmodel.modelParams.NB,len);
        v_RNEA = cell(len,1);
        resRNEA.tau = zeros(length(dataTime),dmodel.NB);
        
        for i = 1 : len
             myRNEA = myRNEA.setState(q(:,i), dq(:,i));
             myRNEA = myRNEA.setY(y_RNEA2(:,i));
             myRNEA = myRNEA.solveID();

             d(:,i) = myRNEA.d;
             
             v_RNEA{i,:} = myRNEA.v; 
             resRNEA.tau(i,:) = myRNEA.tau;
        end
        
        clear fx;
        
        %% Build bias b_Y manually
        % b_Y has to be consistent with Ymatrix
        
        load('./experiments/humanFixedBase/data/subjectSizeParams.mat');

        currentParams = subjectParams(subjectID);

        footMass =  currentParams.footMass;
        posP_0 = [0; 0; (0.5*currentParams.footHeight)];
        footIxx =  currentParams.footIxx;
        footIyy =  currentParams.footIyy;
        footIzz =  currentParams.footIzz;

        b_Y_3sens = zeros (size(data.y_3sens)); 
        R_imu_2 = X_imu_2(1:3,1:3);

        a_G_grav = [0;0;0;0;0;-9.8100]; %Featherstone-like notation, in global reference
        X_G_0 = currentTrialSens.X_G_0;
        X_0_G = InverseAdjTransform(X_G_0);
        a_0_grav = X_0_G * a_G_grav;

        I_0 = createSpatialInertia(footIxx,footIyy,footIzz,footMass,posP_0);

        b_Y_3sens(1:6,1:len)   = repmat((-XStar_fp_0 * I_0 * a_0_grav),1,len);
        
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
            b_Y_3sens(10:12,i) = cross(A,B);
        end 
        
        clear A;
        clear B;
        clear footIxx;
        clear footIyy;
        clear footIzz;
        
        
        
%         % exploiting velocitiy v_RNEA coming from RNEA class computation
%         for i = 1 : len 
%             A =R_imu_2*v_RNEA{i,1}(1:3,2);
%             B =((X_imu_2(4:6,1:3)*v_RNEA{i,1}(1:3,2))+(R_imu_2*v_RNEA{i,1}(4:6,2)));
%             b_Y_3sens(4:6,i) = cross(A,B);
%         end 

        %% Computing MAP method
        fprintf('\nMAP computation with 3 sensors\n');
        resMAP_3sens.d  = zeros(26*dmodel.NB,len);
        resMAP_3sens.Sd = cell(len,1);
        resMAP_3sens.Ymatrix = cell(len,1);
        resMAP_3sens.b_Y = zeros (myMAP_3sens.IDsens.sensorsParams.m,len);
        
        for i = 1 : len

            myMAP_3sens = myMAP_3sens.setState(data.q(i,:)', data.dq(i,:)');
            myMAP_3sens = myMAP_3sens.setY(data.y_3sens(:,i));
            myMAP_3sens = myMAP_3sens.setYmatrix(Y_3sens);
            myMAP_3sens = myMAP_3sens.setBias(b_Y_3sens(:,i));
            myMAP_3sens = myMAP_3sens.solveID();

            resMAP_3sens.d(:,i)       = myMAP_3sens.d;
            resMAP_3sens.Sd{i,1}      = myMAP_3sens.Sd;
            resMAP_3sens.Ymatrix{i,1} = myMAP_3sens.IDsens.sensorsParams.Y; 
            resMAP_3sens.b_Y(:,i)     = myMAP_3sens.IDsens.sensorsParams.bias;
            
             if mod(i-1,100) == 0
                    fprintf('Processing %d %% of the dataset\n', round(i/len*100));
             end
        end
        
        % ========end MAP
        %% Rearrange solution

        for i = 1 : dmodel.NB
             
            link = strrep(myMAP_3sens.IDmodel.modelParams.linkname{i}, '+', '_');
            joint = strrep(myMAP_3sens.IDmodel.modelParams.jointname{i}, '+', '_');
             
            % initialize variables
            eval(['resMAP_3sens.a_'    link ' = zeros(6, len );']);
            eval(['resMAP_3sens.Sa_'    link ' = cell(len, 1);']);
            
            eval(['resMAP_3sens.fB_'    link ' = zeros(6, len );']);
            eval(['resMAP_3sens.SfB_'    link ' = cell(len, 1);']);
            
            eval(['resMAP_3sens.f_'    link ' = zeros(6, len );']);
            eval(['resMAP_3sens.Sf_'     link ' = cell(len, 1);']);
            
            eval(['resMAP_3sens.tau_'    joint ' = zeros(1, len );']);
            eval(['resMAP_3sens.Stau_'    joint ' = zeros(len, 1);']);
            
            eval(['resMAP_3sens.fx_'    link ' = zeros(6, len );']);
            eval(['resMAP_3sens.Sfx_'    link ' = cell(len, 1);']);
            
            eval(['resMAP_3sens.d2q_'    joint ' = zeros(1, len );']);
            eval(['resMAP_3sens.Sd2q_'    joint ' = zeros(len, 1);']);
             
            
            for j = 1 : len
                 %a
                 ind  = '1 + 26*(i-1) : 26*(i-1) +  6';
                 eval(['resMAP_3sens.a_'    link '(:,j) = resMAP_3sens.d      (' ind '       ,j);'])
                 eval(['resMAP_3sens.Sa_'   link '{j,1} = resMAP_3sens.Sd{j,1}(' ind ',' ind ' );'])
                 %fB
                 ind  = '7 + 26*(i-1) : 26*(i-1) + 12';
                 eval(['resMAP_3sens.fB_'   link '(:,j) = resMAP_3sens.d      (' ind '        ,j);'])
                 eval(['resMAP_3sens.SfB_'  link '{j,1} = resMAP_3sens.Sd{j,1}(' ind ',' ind '  );'])
                 %f
                 ind  = '13 + 26*(i-1) : 26*(i-1) + 18';
                 eval(['resMAP_3sens.f_'    link '(:,j) = resMAP_3sens.d      (' ind '        ,j);'])
                 eval(['resMAP_3sens.Sf_'   link '{j,1} = resMAP_3sens.Sd{j,1}(' ind ',' ind '  );'])
                 %tau
                 ind  = '19 + 26*(i-1) : 26*(i-1) + 19';
                 eval(['resMAP_3sens.tau_'  joint '(:,j) = resMAP_3sens.d      (' ind '        ,j);'])
                 eval(['resMAP_3sens.Stau_' joint '(j,1) = resMAP_3sens.Sd{j,1}(' ind ',' ind '  );'])
                 %fx
                 ind  = '20 + 26*(i-1) : 26*(i-1) + 25';
                 eval(['resMAP_3sens.fx_'   link '(:,j) = resMAP_3sens.d      (' ind '        ,j);'])
                 eval(['resMAP_3sens.Sfx_'  link '{j,1} = resMAP_3sens.Sd{j,1}(' ind ',' ind '  );'])
                 %ddq
                 ind  = '26 + 26*(i-1) : 26*(i-1) + 26';
                 eval(['resMAP_3sens.d2q_'  joint '(:,j) = resMAP_3sens.d      (' ind '        ,j);'])
                 eval(['resMAP_3sens.Sd2q_' joint '(j,1) = resMAP_3sens.Sd{j,1}(' ind ',' ind '  );'])   
             end
        end 
      
        
        
        %% =========================== MAP with 4 sensors(ddq+ftx(1&2)+fp+imu) ===============================
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
        data.y_4sens  = [data.y_4sens, data.ddq];
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

        
        %% Build bias b_Y manually
        % b_Y has to be consistent with Ymatrix

        load('./experiments/humanFixedBase/data/subjectSizeParams.mat');

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

            myMAP_4sens = myMAP_4sens.setState(data.q(i,:)', data.dq(i,:)');
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
        
        finalResults(subjectID,trialID).dataTime = dataTime;
        finalResults(subjectID,trialID).resMAP.Sprior_inv = [  zeros(38),          zeros(38,14)       ;
                                                              zeros(14,38), full(dmodel.Sw_inv.matrix)];
        finalResults(subjectID,trialID).resMAP.D = myMAP_2sens.D;
        finalResults(subjectID,trialID).resMAP.b_D = myMAP_2sens.b_D;
        finalResults(subjectID,trialID).resMAP.SD_inv = full(dmodel.Sv_inv.matrix);
                                                          
        finalResults(subjectID,trialID).resMAP_2sens = resMAP_2sens;
        finalResults(subjectID,trialID).resMAP_2sens.m = mySens_2sens.m;
        finalResults(subjectID,trialID).resMAP_2sens.data.y = data.y_2sens;
        finalResults(subjectID,trialID).resMAP_2sens.Sy_inv = full(mySens_2sens.sensorsParams.Sy_inv);
        finalResults(subjectID,trialID).resMAP_2sens.jIndex = myMAP_2sens.jIndex;
        
        finalResults(subjectID,trialID).resMAP_3sens = resMAP_3sens;
        finalResults(subjectID,trialID).resMAP_3sens.m = mySens_3sens.m;
        finalResults(subjectID,trialID).resMAP_3sens.data.y = data.y_3sens;
        finalResults(subjectID,trialID).resMAP_3sens.Sy_inv = full(mySens_3sens.sensorsParams.Sy_inv);
        finalResults(subjectID,trialID).resMAP_2sens.jIndex = myMAP_3sens.jIndex;
        
        finalResults(subjectID,trialID).resMAP_4sens = resMAP_4sens;
        finalResults(subjectID,trialID).resMAP_4sens.m = mySens_4sens.m;
        finalResults(subjectID,trialID).resMAP_4sens.data.y = data.y_4sens;
        finalResults(subjectID,trialID).resMAP_4sens.Sy_inv = full(mySens_4sens.sensorsParams.Sy_inv);
        finalResults(subjectID,trialID).resMAP_2sens.jIndex = myMAP_4sens.jIndex;

%         finalResults(subjectID,trialID).resMAP.data.Sy = data.Sy;

         clear res;

    end
    fprintf('\n');
end
%% storing results
save('./experiments/humanFixedBase/intermediateDataFiles/finalResults.mat','finalResults');

fprintf('---------\n');
fprintf('Done!\n');
