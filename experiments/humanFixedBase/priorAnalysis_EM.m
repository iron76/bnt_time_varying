clear;clc;close all;

onlyIMUMeasurement = 'off';

if(strcmp(onlyIMUMeasurement,'on')==1)
    load('./experiments/humanFixedBase/intermediateDataFiles/savedBNet_onlyIMU.mat');
    fprintf('EM learning on the only IMU dataset\n');
else
    load('./experiments/humanFixedBase/intermediateDataFiles/savedBNet_IMUandFTS.mat');
    fprintf('EM learning on the IMU and FTS dataset\n');
end

fprintf('Loaded the saved Bayesian network \n Nodes to be learned using EM : \n');

NB = myModel.modelParams.NB;
    i_learn = [];
    for i = 1 : length(myBNEA.bnt.nodes.labels)
       if(strcmp(onlyIMUMeasurement,'on')==1) 
           if strcmp(myBNEA.bnt.nodes.labels{i}(1:2), 'y_')...
                   && ~strcmp(myBNEA.bnt.nodes.labels{i}(end-2:end), 'ftx')...
                   &&  ~strcmp(myBNEA.bnt.nodes.labels{i}(end-2:end), 'fts')
              fprintf('%s\n',myBNEA.bnt.nodes.labels{i}); 
              i_learn = [i_learn, i];
           end
       else
           if strcmp(myBNEA.bnt.nodes.labels{i}(1:2), 'y_')...
                   && ~strcmp(myBNEA.bnt.nodes.labels{i}(end-2:end), 'ftx')%...
              %     &&  ~strcmp(myBNEA.bnt.nodes.labels{i}(end-2:end), 'fts')
              fprintf('%s\n',myBNEA.bnt.nodes.labels{i});
              i_learn = [i_learn, i];
           end
       end
    end
    
    figure;
    

    fprintf('Learning EM Bnet parameters\n');
    % This parameters should be kept small (1e-15)for
    % simulated data. For real data is plays the role
    % of a regularization factor and should be tuned so
    % as to have non-decresing EM steps.
    cov_prior_weight = cell(18,1);
    for i = 1:length(i_learn);
        cov_prior_weight{i_learn(i)} = 1e-3;%1e1 ;
    end
    
   % standardisation matrix for IMU in the 0 frame 
   SimuAng_2 = [1e0 1e-3 1e0];
   SimuLin_2 = [1e-3 1e0 1e-3];
    
   % standardisation matrix for FT in the 0 frame
   SfTStandMu_0 = [1e0 1e0 1e-3];
   SfTStandF_0 = [1e-3 1e-3 1e0];
    
    %SfTStandMu_0 = [1e-10 1e-3 1e-10];
    %SfTStandF_0 = [1e-3 1e-10 1e-3];
  
   cov_prior_weight{17} = diag([SfTStandMu_0,SfTStandF_0]);
   cov_prior_weight{9} = diag([SimuAng_2 SimuLin_2]); 
   
   %% standardisation applied in FT frame
 %   cov_prior_weight{17} = (XStar_fp_0^(-1))*diag([SfTStandMu_0 SfTStandMu_0]);
    maxSteps = 6;
    
    [bnetHat, ll] = EM_bnet_learn(bnet, sample, cov_prior_weight, i_learn, maxSteps);

    dir_ind = cell2mat(myBNEA.bnt.nodes.index);
    inv_ind(dir_ind) = 1:length(dir_ind);
    
    fprintf('Analysing results\n');

    % Depending on the number of samples 'n' the updates can
    % be quite relevant. With increasing 'n' the magnitude of
    % the updates tends to become smaller.


    fprintf('Analysing and plotting results \n');
    
if(strcmp(onlyIMUMeasurement,'on')==1)
    save('./experiments/humanFixedBase/intermediateDataFiles/EMResult_onlyIMU.mat','myModel','mySens','bnet','bnetHat','ll','i_learn');
    
else
    save('./experiments/humanFixedBase/intermediateDataFiles/EMResult_IMUandFTS.mat','myModel','mySens','bnet','bnetHat','ll','i_learn');
end

priorAnalysis_plotResults(onlyIMUMeasurement,myBNEA,ll,bnetHat,i_learn,mySens);

    