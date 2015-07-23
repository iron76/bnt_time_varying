function res = testLearnBNEAIP(dmodel, ymodel)

res = 0;

n  = 10;         %number of samples
NB = dmodel.NB;  %number of links

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myBNEAIP  = BNEAIP(myModel, mySens);
myBNEAIP  = myBNEAIP.setEngine('jtree_inf_engine');


bnet    = cell(1,n);
engine  = cell(1,n);
sample  = cell(NB*6+myBNEAIP.IDsens.sensorsParams.ny,n);

for i = 1 : n
   q  = rand(dmodel.NB,1);
   dq = rand(dmodel.NB,1);
   
   myBNEAIP = myBNEAIP.setState(q, dq);
   [y,yc] = myBNEAIP.simY();
   myBNEAIP = myBNEAIP.setY(y);
   
   bnet{1,i}   = myBNEAIP.bnt.bnet;
   engine{1,i} = myBNEAIP.bnt.engine;
   sample(:,i) = yc;
end

% This parameters should be kept small (1e-15)for
% simulated data. For real data is plays the role
% of a regularization factor and should be tuned so
% as to have non-decresing EM steps.

cov_prior_weight = 1e-15;
[bnetHat, ll] = EM_bnet_learn(bnet, sample, cov_prior_weight);

if (sum(diff(ll) < 0) ~= 0) && (norm(ll(end-1)-ll(end)) > 1e-3)
   disp('Something wrong with the EM algorithm. Declaring the test failed!')
   res = 1;
end

dir_ind = cell2mat(myBNEAIP.bnt.nodes.index);
inv_ind(dir_ind) = 1:length(dir_ind);

for i = 1 : length(dir_ind(1:NB*6))
   cov_ini = get_field(myBNEAIP.bnt.bnet.CPD{dir_ind(i)}, 'cov');
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
   cov_ini = get_field(myBNEAIP.bnt.bnet.CPD{dir_ind(i)}, 'cov');
   cov_est = get_field(bnetHat.CPD{dir_ind(i)},         'cov');
   cov_upd = cov_est - cov_ini;
   fprintf('[INFO] %s was updated by %f \n', myBNEAIP.bnt.nodes.labels{i}, norm(cov_upd)./norm(cov_ini));
end
