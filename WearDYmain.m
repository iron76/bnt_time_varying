
clear;clc;close all;

% forcing the casual generator to be const
rng(1);

%% add\create folders

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

%% add\create files .mat

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
subjectList = 1;
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

        %% build model valid for all MAP case 
        load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');

        currentModel = humanThreeLinkModelFromURDF(subjectID);

        humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
        humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 
        
        dmodel  = currentModel.dmodel;                      % deterministic model
        dmodel  = autoTreeStochastic(dmodel, 1e-4);         % probabilistic model for D equation (added Sv and Sw)
        myModel = model(dmodel);    
        
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
         %% Computing v
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
        
        %% Computing EM
        fprintf('\nStarting EM procedure\n');    
        
        sigma_ygivend_init = inv(full(myMAP_4sens.IDsens.sensorsParams.Sy_inv));
        
        [sigma_ygivend,Sigma_EM] = computeEM_test(myMAP_4sens,data.q',data.dq' ,Ymatrix_4sens, data.y_4sens,b_Y_4sens, sigma_ygivend_init, 'MAX_ITER', 4);

        % ========end EM
    end
    fprintf('\n');
end
%% storing results
%save('./experiments/humanFixedBase/intermediateDataFiles/finalResults.mat','finalResults');

fprintf('---------\n');
fprintf('Done!\n');


%% MAP