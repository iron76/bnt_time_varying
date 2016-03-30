clear;
close all;

%% selected subjects and trials
subjectList = 1;
trialList = 1;  

dispSensors = 'off'; % turn on to see sensor data in link and in sensor frames
simpleMeasurement = 'off'; % use only IMU
selectedPercentage = 5; % percentage of time points used for estimation (reduce to speed up while losing accuracy)

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
    
    
    if(strcmp(simpleMeasurement,'on')==1)
        %% =====structure from files
        data.parts    = {'torso'};
        data.labels   = {'imu'  };
        data.ndof     = {6};
        data.index    = {'1:6'};

        %% =====structure of sensors
        sens.parts    = {'torso'};       %force of the forceplate is trasmitted into the leg
        sens.labels   = {'imu'};
        sens.ndof     = {6};
        label_to_plot = {'imu'};
    else
                %% =====structure from files
        data.parts    = {'leg'         ,'torso'};
        data.labels   = {'fts'         ,'imu'  };
        data.ndof     = {6             ,6      };
        data.index    = {'1:6'         ,'1:6'  };

        %% =====structure of sensors
        sens.parts    = {'leg'         ,'torso'};       %force of the forceplate is trasmitted into the leg
        sens.labels   = {'fts'         ,'imu'  };
        sens.ndof     = {6             ,6      };

        label_to_plot = {'fts'         ,'imu'  };
    end
  
    %% build models
    load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');
    
    currentModel = humanThreeLinkModelFromURDF(subjectID);
    
    humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
    humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 
   
    dmodel  = currentModel.dmodel;                      %deterministic model 
    ymodel  = humanThreeLinkSens(dmodel, sens);  
   
    dmodel  = autoTreeStochastic(dmodel, 1e-6);   % probabilistic model for D equation (added Sv and Sw)
    ymodel  = humanThreeLinkSensStochastic(ymodel);     % probabilistic model for Y(q,dq) d = y (added Sy)
   
    myModel = model(dmodel);
    mySens  = sensors(ymodel);  

%     myMAP   = MAP(myModel, mySens);
    %% ================================ RNEA ==================================
    % Computing d using Newton-Euler with inverse dynamics.
    % Function ID was modified for having velocities.  So the path is not:
    % ../bnt_time_varying/extern/featherstone/dynamics/ID.m
    % but
    % ../bnt_time_varying/experiments/humanFixedBase/helperFunction/IDv.m

    tau = zeros(size(q))';
    a = cell (size(q))';
    fB = cell (size(q))';
    f = cell (size(q))';
    fx = zeros (6,1);

    v = cell(size(q))';

    d_temp = zeros(26*dmodel.NB,1);
    d = zeros (26*dmodel.NB, len);
    
    fext    = cell(1,2);
    for i = 1 : dmodel.NB
         fext{i}    = fx;
    end

    for i = 1:len
    
        [tau_i, a_i,v_i,fB_i, f_i] = IDv(dmodel, q(:,i), dq(:,i), ddq(:,i), fext);
        tau(i,:) = tau_i;
        a(i,:) = a_i;
        fB(i,:) = fB_i;
        f(i,:) = f_i;  
        v(i,:) = v_i;
      
        for j = 1 : dmodel.NB
             d_temp((1:26)+(j-1)*26) = [a_i{j}; fB_i{j}; f_i{j}; tau(i,j); fx; ddq(j,i)];
        end
      
        d(:,i) = d_temp;
    end

    fRNEA = cell2mat(f(:,1)');
    aRNEA = cell2mat(a(:,2)');
    ddqRNEA = ddq;
    
    clear d_temp;
    clear tau_i; 
    clear a_i;
    clear fB_i;
    clear f_i;
    clear v_i;
    
    %% Build data.y anda data.Sy 
    % data.y are ordered in:  
    % - angular-linear notation 
    % - the form [f1 a2 ftx1 ftx2 ddq1 ddq2]
    % - sensor frame

    %% ==== data.y (in sensor frame)
    data.y  = [];
    for i = 1 : length(sens.labels)
         eval(['data.y  = [data.y  data.ys_sensFrame_' sens.labels{i} '];']);
    end
    %% Add null external forces ftx = 0
    data.y  = [data.y, zeros(len,6*dmodel.NB)];
    %% Add ddq measurements
    data.y  = [data.y, data.ddq];
    data.y = data.y';

    %% ==== data.Sy
    data.Sy = [];
    
     
   
%     for i = 1 : length(myMAP.IDsens.sensorsParams.labels)
% 
%          data.Sy = [data.Sy; diag(myMAP.IDsens.sensorsParams.Sy{i})];
%          
%          if(isempty(strfind(myMAP.IDsens.sensorsParams.labels{i},'d2q')))
%              figure(ceil(i/4));
%              subplot(2,2,mod(i,4)+1);
%              bar(myMAP.IDsens.sensorsParams.Sy{i},'stacked');
%              title(myMAP.IDsens.sensorsParams.labels{i});
%          end
%     end
%     data.Sy = repmat(data.Sy, 1, len-1);
%     data.Sy = [data.Sy data.Sy(:,end)];
%     drawnow();

%     %% ================================ MAP ===================================
%     %% Build Ymatrix manually
%     % Ymatrix has to be consistent with measurements form [f1 a2 ftx1 ftx2 ddq1 ddq2]
%     
%     Y = zeros (ymodel.m,26*dmodel.NB);
%     Y(10:12,27:32) = X_imu_2(4:6,:);
%     Y(13:18,20:25) = eye(6);
%     Y(19:24,46:51) = eye(6);
%     Y(25,26) = eye(1);
%     Y(26,52) = eye(1);
% 
%     
%     Ymatrix = cell(len,1);
%     for i = 1 : len
%         %the only row in Ymatrix that is time varying
%         Y(1:6,13:18) = XStar_fp_0 * XStar_0_1{i};
%         Ymatrix{i} = Y; 
%     end
%     
%     clear Y;
%    
    %% Build bias b_Y manually
    % b_Y has to be consistent with Ymatrix
    
    load('./experiments/humanFixedBase/data/subjectSizeParams.mat');
    
    currentParams = subjectParams(subjectID);
    
    footMass =  currentParams.footMass;
    posP_0 = [0; 0; (0.5*currentParams.footHeight)];
    footIxx =  currentParams.footIxx;
    footIyy =  currentParams.footIyy;
    footIzz =  currentParams.footIzz;
    
    b_Y = zeros (size(data.y)); 
    R_imu_2 = X_imu_2(1:3,1:3);
    
    a_G_grav = [0;0;0;0;0;-9.8100]; %Featherstone-like notation, in global reference
    X_0_G = currentTrialSens.X_G_0';
    a_0_grav = X_0_G * a_G_grav;
          
    I_0 = createSpatialInertia(footIxx,footIyy,footIzz,footMass,posP_0);

    b_Y(1:6,1:len)   = repmat((-XStar_fp_0 * I_0 * a_0_grav),1,len);
    
    for i = 1 : len   
        A =R_imu_2*v{i,2}(1:3,1);
        B =((X_imu_2(4:6,1:3)*v{i,2}(1:3,1))+(R_imu_2*v{i,2}(4:6,1)));
        b_Y(10:12,i) = cross(A,B);
    end 

    clear A;
    clear B;
    
        %% displaying sensorCovariances
    for i = 1 : length(mySens.sensorsParams.labels)

         %data.Sy = [data.Sy; diag(mySens.sensorsParams.Sy{i})];
         if(isempty(strfind(mySens.sensorsParams.labels{i},'d2q')))
             figure(ceil(i/4));
             subplot(2,2,mod(i,4)+1);
             imagesc(mySens.sensorsParams.Sy{i});
             title([mySens.sensorsParams.labels{i},' sensorFrame']);
         end
    end
    
    
    %% recomputing data in link frames
    
    if(strcmp(simpleMeasurement,'on')==1)
        data.yBiasComp = data.y;
        data.yLinkFrame = data.yBiasComp;
        data.yLinkFrame(1:6,:) = ((X_imu_2)^(-1))*data.yBiasComp(1:6,:);
        yLinkFrameRNEA = [ aRNEA;...
            zeros(12,length(data.dataTime));ddqRNEA];
    else
        data.yBiasComp = zeros(size(data.y));
        data.yLinkFrame = zeros(size(data.y));
        yLinkFrameRNEA = zeros(size(data.y));

        tTot = length(data.dataTime);
        for i = 1:tTot
            %% data must be brought to the link-wise frames
           %% first we combine bias matrix (in sensor frame) with sensor measurements

           data.yBiasComp(:,i) = data.y(:,i);% - b_Y(:,i);

           %% recomputing the Force plate and IMU data into link frames

           data.yLinkFrame(1:6,i) = ((XStar_fp_0 * XStar_0_1{i} ))^(-1)*data.yBiasComp(1:6,i);
           data.yLinkFrame(7:12,i) = ((X_imu_2)^(-1))*data.yBiasComp(7:12,i);
           data.yLinkFrame(13:end,i) = data.yBiasComp(13:end,i);
           
        end
        %% rotating the senor covariances to bring to link frame (only for FTS and IMU) and only in initial configuration
      %% future analysis
           
           mySens.sensorsParams.sensXLink{1} = XStar_fp_0 * XStar_0_1{1};
           mySens.sensorsParams.sensXLink{2} = X_imu_2;
           mySens.sensorsParams.Sy{1} = ((XStar_fp_0 * XStar_0_1{1} ))^(-1)*mySens.sensorsParams.Sy{1}*(((XStar_fp_0 * XStar_0_1{1} ))^(-1))';
           mySens.sensorsParams.Sy{2} =  ((X_imu_2)^(-1))*mySens.sensorsParams.Sy{2}*((X_imu_2)^(-1))';
           
       my = 1;idSy_inv = []; jdSy_inv = []; dSy_inv=[];    
       %% recomputing sensorParams.sy_inv    
        for i = 1 : mySens.sensorsParams.ny
            dy = mySens.sensorsParams.sizes{i,1};
           [ii, jj, ss] = placeSubmatrixSparse(my, my, inv(mySens.sensorsParams.Sy{i,1}));
           idSy_inv = [idSy_inv; ii];
           jdSy_inv = [jdSy_inv; jj];
           dSy_inv  = [dSy_inv;  ss];
           my = my + dy;
        end
        mySens.sensorsParams.Sy_inv = sparse(idSy_inv, jdSy_inv, dSy_inv);
           
          % mySens.sensorsParams.Sy_inv{1} = mySens.sensorsParams.Sy{1}^(-1);
          % mySens.sensorsParams.Sy_inv{2} = mySens.sensorsParams.Sy{2}^(-1);
    
        yLinkFrameRNEA = [ fRNEA;aRNEA;...
            zeros(12,length(data.dataTime));ddqRNEA];
    end
    
    
    %% displaying sensorCovariances
    for i = 1 : length(mySens.sensorsParams.labels)

         %data.Sy = [data.Sy; diag(mySens.sensorsParams.Sy{i})];
         if(isempty(strfind(mySens.sensorsParams.labels{i},'d2q')))
             figure(1+ceil(i/4));
             subplot(2,2,mod(i,4)+1);
             imagesc(mySens.sensorsParams.Sy{i});
             title([mySens.sensorsParams.labels{i},'linkFrame']);
         end
    end
 
    
    %% Creating the Bayesian Network

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    myBNEA  = BNEA(myModel, mySens);
    myBNEA  = myBNEA.setEngine('jtree_inf_engine');
    NB      = myModel.modelParams.NB;
    [~, nTot]  = size(data.y);

    
    n = round(nTot*(selectedPercentage/100));
    
    bnet    = cell(1,n);
    engine  = cell(1,n);
    sample  = cell(NB*6+myBNEA.IDsens.sensorsParams.ny,n);
    sampleTest  = cell(NB*6+myBNEA.IDsens.sensorsParams.ny,n);
    
    printPts = floor(linspace(1,n,10));
    printPts = [printPts(2:end) printPts(end)];
    printPtsCtr = 1;
    fprintf('Updating Bayesian network with state and Y using %d points \n',n);
    %disp(bnet);
    
    sortYidxComplexMeasurement = [4:6 1:3 10:12 7:9 13:26]; %% hack to fix evidence ordering
    %sortYidxComplexMeasurement = 1:26;%[4:6 1:3 10:12 7:9 13:26];
    sortYidxSimpleMeasurement = [4:6 1:3 7:20]; 
    
    for i = 1 : n

      if(i == printPts(printPtsCtr))
          fprintf('Completed %d%%\n',printPtsCtr*10);
          printPtsCtr = printPtsCtr+1;
      end
      
       myBNEA = myBNEA.setState(data.q(i,:)', data.dq(i,:)');
       %myBNEA = myBNEA.setY(data.yLinkFrame(:,i));
       
       
       if(strcmp(simpleMeasurement,'on')==1)
          myBNEA=myBNEA.setY(data.yLinkFrame(sortYidxSimpleMeasurement,i));  
       else
          myBNEA=myBNEA.setY(data.yLinkFrame(sortYidxComplexMeasurement,i));
       end   
       %myBNEA = myBNEA.setY(yLinkFrameRNEA(:,i)); %% set from RNEA
       % myBNEA = myBNEA.setY(zeros(size(data.yLinkFrame(:,i)))); %% set
       % all zeros

       bnet{1,i}   = myBNEA.bnt.bnet;
       engine{1,i} = myBNEA.bnt.engine;
       sample(:,i) = myBNEA.evd';
        
       %% sampleTest is used to check accuracy of stored evidence in sample
       if(strcmp(simpleMeasurement,'on')==1)
           sampleTest(2,i) = {data.yLinkFrame(20,i)}; % ddq2
           sampleTest(4,i) = {data.yLinkFrame(13:18,i)}; % fx2
           sampleTest(6,i) = {data.yLinkFrame(19,i)}; % ddq1
           sampleTest(9,i) = {data.yLinkFrame(1:6,i)}; % aimu
           sampleTest(15,i) = {data.yLinkFrame(7:12,i)}; % fx1
       else
           %% sampleTest
           sampleTest(2,i) = {data.yLinkFrame(26,i)}; % ddq2
           sampleTest(4,i) = {data.yLinkFrame(19:24,i)}; % fx2
           sampleTest(6,i) = {data.yLinkFrame(25,i)}; % ddq1
           sampleTest(9,i) = {data.yLinkFrame(7:12,i)}; % aimu
           sampleTest(15,i) = {data.yLinkFrame(13:18,i)}; % fx1
           sampleTest(17,i) = {data.yLinkFrame(1:6,i)}; % fts
       end
    end
    end
end

if(strcmp(dispSensors,'on')==1)
 %   plotMeasurements(data.dataTime,data.y,'sensorFrame',...
  %                                 data.yLinkFrame,'LinkFrame');
   % plotMeasurements(data.dataTime(1:n),data.yLinkFrame(:,1:n),'sensorFrame',...
   %                                cell2mat(sample),'evidence');
  % plotMeasurements(data.dataTime,data.yLinkFrame,'SensorLinkFrame',...
  %                                   yLinkFrameRNEA,'RNEALinkFrame');
   plotEvidence(data.dataTime(1:n),cell2mat(sample),cell2mat(sampleTest));                               
                               
end
save('./experiments/humanFixedBase/intermediateDataFiles/savedBNet.mat','data','myModel','mySens','myBNEA','bnet','engine','sample','data','XStar_fp_0');