
% checkRNEA_ID
% Script to compute vector d = [d_1,d_2,...,d_NB] comparing two
% equivalent methods: 
% - METHOD 1: using Featherstone ID --> d
% - METHOD 2: using RNEA class --> d_RNEA

clear; close all; clc;

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
    
    
    % =====structure from files
    data.parts    = {'leg'         ,'torso'};
    data.labels   = {'fts'         ,'imu'  };
    data.ndof     = {6             ,6      };
    data.index    = {'1:6'         ,'1:6'  };

    % =====structure of sensors
    sens.parts    = {'leg'         ,'torso'};       %force of the forceplate is trasmitted into the leg
    sens.labels   = {'fts'         ,'imu'  };
    sens.ndof     = {6             ,6      };

    label_to_plot = {'fts'         ,'imu'  };
    

    %% build models
    load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');
    
    currentModel = humanThreeLinkModelFromURDF(subjectID);
    
    humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
    humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 
   
    dmodel  = currentModel.dmodel;                      %deterministic model 
    ymodel  = humanThreeLinkSens(dmodel, sens);  
   
    dmodel  = autoTreeStochastic(dmodel, 1e-5, 1e-4);   % probabilistic model for D equation (added Sv and Sw)
    ymodel  = humanThreeLinkSensStochastic(ymodel);     % probabilistic model for Y(q,dq) d = y (added Sy)
   
    myModel = model(dmodel);
    mySens  = sensors(ymodel); 

    %% ======METHOD 1: Computing d using Newton-Euler with Featherstone ID
    tic;

    tau = zeros(size(q))';
    a = cell (size(q))';
    fB = cell (size(q))';
    f = cell (size(q))';
    fx = zeros (6,1);

    d_temp = zeros(26*dmodel.NB,1);
    d = zeros (26*dmodel.NB, len);

    fext    = cell(1,2);
    for i = 1 : dmodel.NB
            fext{i}    = fx;
    end

    for i = 1:len
    
        [tau_i, a_i,fB_i, f_i] = ID(dmodel, q(:,i), dq(:,i), ddq(:,i), fext);
        tau(i,:) = tau_i;
        a(i,:) = a_i;
        fB(i,:) = fB_i;
        f(i,:) = f_i;  
      
        for j = 1 : dmodel.NB
             d_temp((1:26)+(j-1)*26) = [a_i{j}; fB_i{j}; f_i{j}; tau(i,j); fx; ddq(j,i)];
        end
      
        d(:,i) = d_temp;
    end

    t_ID = toc;
    disp(['CPU time for d computation with ID method is: ' num2str(t_ID) '[sec]']);
    disp('/\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ ')
 
    clear d_temp;
    clear tau_i; 
    clear a_i;
    clear fB_i;
    clear f_i;

    %% =====METHOD 2: Computing d using Newton-Euler with RNEA method 
    tic;

    ymodel_RNEA  = autoSensRNEA(dmodel);
    mySens_RNEA  = sensors(ymodel_RNEA);
    myRNEA       = RNEA(myModel, mySens_RNEA);
      
    y_RNEA_f = zeros(6*dmodel.NB, len);
    y_RNEA_ddq = zeros(dmodel.NB, len);
    fx = cell(dmodel.NB);
  
    %Ordering y_RNEA in the form [fx1 fx2 ddq1 ddq2]
    for i = 1 : dmodel.NB
        for t = 1 : len
            fx{i,1} = zeros(6,1); 
            y_RNEA_f(((1:6)+(i-1)*6), t) = [fx{i,1}];
            y_RNEA_ddq(i, t) = [ddq(i,t)];
        end
        y_RNEA = [y_RNEA_f ; y_RNEA_ddq];
    end


    d_RNEA = zeros (26*myRNEA.IDmodel.modelParams.NB,len);
        for i = 1 : len
             myRNEA = myRNEA.setState(q(:,i), dq(:,i));
             myRNEA = myRNEA.setY(y_RNEA(:,i));
             myRNEA = myRNEA.solveID();
       
             d_RNEA(:,i) = myRNEA.d; 
        end
   
    t_RNEA = toc;
    disp(['[1st] CPU time for d computation with RNEA class method is: ' num2str(t_RNEA) '[sec]']);

    %% =====check

    if (sum(d-d_RNEA) ~= 0);
        disp('Something wrong with d computation. Check methods.');
        res = 1; 
    else
        disp('Methods 1 and 2 are equivalent.')
    end

   end
   fprintf('\n');
end