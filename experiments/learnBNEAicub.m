clear all
close all
load preprocess2.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myBNEA  = BNEA(myModel, mySens);
myBNEA  = myBNEA.setEngine('jtree_inf_engine');
NB      = myModel.modelParams.NB;
[~, n]  = size(y);

bnet    = cell(1,n);
engine  = cell(1,n);
sample  = cell(NB*6+myBNEA.IDsens.sensorsParams.ny,n);

for i = 1 : n
  
   myBNEA = myBNEA.setState(data.q(:,i), data.dq(:,i));
   myBNEA = myBNEA.setY(y(:,i));
   
   bnet{1,i}   = myBNEA.bnt.bnet;
   engine{1,i} = myBNEA.bnt.engine;
   sample(:,i) = myBNEA.evd';
end

i_learn = [];
for i = 1 : length(myBNEA.bnt.nodes.labels)
   if strcmp(myBNEA.bnt.nodes.labels{i}(1:2), 'y_') && ~strcmp(myBNEA.bnt.nodes.labels{i}(end-2:end), 'ftx')
      i_learn = [i_learn, i];
   end
end

% This parameters should be kept small (1e-15)for
% simulated data. For real data is plays the role
% of a regularization factor and should be tuned so
% as to have non-decresing EM steps.
cov_prior_weight = 0.05;
[bnetHat, ll] = EM_bnet_learn(bnet, sample, cov_prior_weight, i_learn);

dir_ind = cell2mat(myBNEA.bnt.nodes.index);
inv_ind(dir_ind) = 1:length(dir_ind);

for i = 1 : length(dir_ind(1:NB*6))
   cov_ini = get_field(myBNEA.bnt.bnet.CPD{dir_ind(i)}, 'cov');
   cov_est = get_field(bnetHat.CPD{dir_ind(i)},         'cov');
   cov_upd = cov_est - cov_ini;
   if cov_upd
      fprintf('[ERROR] Something wrong with clamped covariance! \n')
      return;
   end
end

% Depending on the number of samples 'n' the updates can
% be quite relevant. With increasing 'n' the magnitude of
% the updates tends to become smaller.

for i = length(dir_ind(1:NB*6))+1 : length(dir_ind)
   cov_ini = get_field(myBNEA.bnt.bnet.CPD{dir_ind(i)}, 'cov');
   cov_est = get_field(bnetHat.CPD{dir_ind(i)},         'cov');
   cov_upd = cov_est - cov_ini;
   fprintf('[INFO] %s was updated by %f \n', myBNEA.bnt.nodes.labels{dir_ind(i)}, norm(cov_upd)./norm(cov_ini));
end

save learn.mat

for i = length(dir_ind(1:NB*6))+1 : length(dir_ind)
   learn.cov_ini{i-length(dir_ind(1:NB*6))} = get_field(myBNEA.bnt.bnet.CPD{dir_ind(i)}, 'cov');
   learn.cov_est{i-length(dir_ind(1:NB*6))} = get_field(bnetHat.CPD{dir_ind(i)},         'cov');
end

save learn_results.mat learn
