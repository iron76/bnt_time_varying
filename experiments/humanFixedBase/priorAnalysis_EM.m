clear;clc;close all;

load('./experiments/humanFixedBase/intermediateDataFiles/savedBNet.mat');
    
fprintf('Loading the saved Bayesian network \n');

NB = myModel.modelParams.NB;
    i_learn = [];
    for i = 1 : length(myBNEA.bnt.nodes.labels)
       if strcmp(myBNEA.bnt.nodes.labels{i}(1:2), 'y_')...
               && ~strcmp(myBNEA.bnt.nodes.labels{i}(end-2:end), 'ftx')%...
          %     &&  ~strcmp(myBNEA.bnt.nodes.labels{i}(end-2:end), 'fts')
          i_learn = [i_learn, i];
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
        cov_prior_weight{i_learn(i)} = 1e-5;%1e1 ;
    end
    
    %cov_prior_weight{3} = 1e-3*eye(6);
    %cov_prior_weight{3}(2,2)= 9e-4;
    %cov_prior_weight{3}(3,3)= 9e-4;
    
    
   SimuAng_2 = [1e0 1e-3 1e0];
   SimuLin_2 = [1e-3 1e0 1e-3];
    
    %% standardisation matrix for FT in the 0 frame
   SfTStandMu_0 = [1e0 1e0 1e-3];
   SfTStandF_0 = [1e-3 1e-3 1e0];
    
    %SfTStandMu_0 = [1e-10 1e-3 1e-10];
    %SfTStandF_0 = [1e-3 1e-10 1e-3];
  % cov_prior_weight{17} = diag([SfTStandMu_0,SfTStandF_0]);
  % cov_prior_weight{9} = diag([SimuAng_2 SimuLin_2]); 
   
    %% standardisation applied in FT frame
 %   cov_prior_weight{17} = (XStar_fp_0^(-1))*diag([SfTStandMu_0 SfTStandMu_0]);
    maxSteps =3;
    
    [bnetHat, ll] = EM_bnet_learn(bnet, sample, cov_prior_weight, i_learn, maxSteps);

    dir_ind = cell2mat(myBNEA.bnt.nodes.index);
    inv_ind(dir_ind) = 1:length(dir_ind);

%     for i = 1:length(ll)
%         for j = ilearn
%            cov_ini = get_field(myBNEA.bnt.bnet.CPD{j}, 'cov');
%            cov_est = get_field(bnetHat.CPD{j},         'cov');
%            cov_upd = cov_est - cov_ini;
%            if cov_upd
%               fprintf('[ERROR] Something wrong with clamped covariance! \n')
%               return;
%            end
%         end
%     end
     
    fprintf('Analysing results\n');

    % Depending on the number of samples 'n' the updates can
    % be quite relevant. With increasing 'n' the magnitude of
    % the updates tends to become smaller.
% 
%     cov_ini = cell(length(i_learn),1);
%     dispVect = i_learn;
%     for j = 1:length(dispVect)
%         cov_ini{j} = get_field(myBNEA.bnt.bnet.CPD{dispVect(j)}, 'cov');
%         cov_ini_diag{j} = diag(get_field(myBNEA.bnt.bnet.CPD{dispVect(j)}, 'cov'));
%     end
%     
%     fprintf('Initial y_fts covariance \n');
%     disp(cov_ini{j});
%     fprintf('EM computed y_fts covariance \n');
%     
%     for j = 1:length(dispVect)
%         s = size(cov_ini{j});
%         cov_est_j = cov_ini{j};
%         cov_est_j_diag = cov_ini_diag{j};
%         
%         for i = 1:length(ll)
%            loc_cov = get_field(bnetHat(i).CPD{dispVect(j)},'cov');
%            loc_cov_diag = diag(loc_cov);
%            cov_est_j = [cov_est_j loc_cov];
%            cov_est_j_diag = [cov_est_j_diag loc_cov_diag];
%            
%            %% printing FT cov
%            if(dispVect(j)==18)
%                fprintf('Step %d, FT Cov diag : \n',i);
%                disp(loc_cov_diag');
%            end
%         end
%         cov_est{j} = cov_est_j;
%         cov_est_diag{j} = cov_est_j_diag;
%    end
%     
%    for i = 1:length(dispVect)
%         figure(i);
%         fprintf('Node %d, cov dim (%d, %d)\n',i,size(cov_est{i}));
%         bar3(1:1+length(ll),cov_est_diag{i}');
%         title(myBNEA.bnt.nodes.labels{i});
%    end

    fprintf('Analysing and plotting results \n');
    save('./experiments/humanFixedBase/intermediateDataFiles/EMResult.mat','myModel','mySens','bnet','bnetHat','ll','i_learn');
    priorAnalysis_plotResults(myBNEA,ll,bnetHat,i_learn,mySens);
    