function [identifiableError,nonIdentifiableError,initialIdentifiableError,initialNotIdentifiableError] = experimentLearnBNEAIP(dmodel, ymodel, seed, nrOfSamples, nrOfIterations)

rng(seed)

n  = nrOfSamples;         %number of samples
NB = dmodel.NB;  %number of links

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the "real model", used to generate the training data
realBNEAIP  = BNEAIP(myModel, mySens);
realBNEAIP  = realBNEAIP.setEngine('jtree_inf_engine');

%disp('Computing identifiable subspace');
[identifiableSubspaceBasis,notIdentifiableSubspaceBasis] = getBaseWrenchIdentifiableSubspace(dmodel);

fprintf('Total mass of original model : %f \n',getTotalMass(myModel.modelParams))

% this is the "training model", on which the EM is run. 
% To check that the inertia parameters are estimated, we see that 
% the total estimated mass is converging on the real values as the EM
% goes on
trainingModel = model(augmentInertialParam(myModel.modelParams,1));

fprintf('Total mass of initial model for training : %f \n',getTotalMass(trainingModel.modelParams))
baseError = identifiableSubspaceBasis'*(inertialParametersFromModel(myModel.modelParams)-inertialParametersFromModel(trainingModel.modelParams));
notIdentifiableError = notIdentifiableSubspaceBasis'*(inertialParametersFromModel(myModel.modelParams)-inertialParametersFromModel(trainingModel.modelParams));

fprintf('Norm of the initial error in base parameters : %f\n', norm(baseError))
fprintf('Norm of the initial error in not identifiable space : %f\n', norm(notIdentifiableError))
fprintf('Norm of the initial error in inertial space (1) : %f\n', norm(inertialParametersFromModel(myModel.modelParams)-inertialParametersFromModel(trainingModel.modelParams)))
fprintf('Norm of the initial error in inertial space (2) : %f\n', norm([baseError;notIdentifiableError]))

initialIdentifiableError = norm(baseError);
initialNotIdentifiableError = norm(notIdentifiableError);


trainingBNEAIP = BNEAIP(trainingModel, mySens);
trainingBNEAIP = trainingBNEAIP.setEngine('jtree_inf_engine');

bnet    = cell(1,n);
engine  = cell(1,n);
sample  = cell(NB*6+realBNEAIP.IDsens.sensorsParams.ny,n);

% generate training data from the real model 
for i = 1 : n
   q  = rand(dmodel.NB,1);
   dq = rand(dmodel.NB,1);
   
   realBNEAIP = realBNEAIP.setState(q, dq);
   [y,yc] = realBNEAIP.simY();
   realBNEAIP = realBNEAIP.setY(y);
   
   trainingBNEAIP = trainingBNEAIP.setState(q, dq);
   
   bnet{1,i}   = trainingBNEAIP.bnt.bnet;
   engine{1,i} = trainingBNEAIP.bnt.engine;
   sample(:,i) = yc;
end

% This parameters should be kept small (1e-15)for
% simulated data. For real data is plays the role
% of a regularization factor and should be tuned so
% as to have non-decresing EM steps.

cov_prior_weight = 1e-15;
[bnetHat, ll] = EM_bnet_learn_IP(bnet, sample, cov_prior_weight, nrOfIterations);

if (sum(diff(ll) < 0) ~= 0) && (norm(ll(end-1)-ll(end)) > 1e-3)
   disp('Something wrong with the EM algorithm. Declaring the test failed!')
   res = 1;
end

dir_ind = cell2mat(realBNEAIP.bnt.nodes.index);
inv_ind(dir_ind) = 1:length(dir_ind);
total_estimated_mass = 0.0;

finalInertialParams = zeros(10*NB,1);
linkToAdd = 1;
for i = 1 : length(dir_ind(1:NB*6))
   cov_ini = get_field(realBNEAIP.bnt.bnet.CPD{dir_ind(i)}, 'cov');
   cov_est = get_field(bnetHat.CPD{dir_ind(i)},         'cov');
   cov_upd = cov_est - cov_ini;
   if cov_upd
      fprintf('[ERROR] Something wrong with clamped covariance! \n')
      return;
   end
   
   % Get total estimated mass
   if( isa(bnetHat.CPD{dir_ind(i)},'gaussian_net_wrench_CPD') )
       inParams = get_field(bnetHat.CPD{dir_ind(i)}, 'inertial_params');
       total_estimated_mass = total_estimated_mass + inParams(1);
       finalInertialParams((1+(linkToAdd-1)*10):linkToAdd*10) = inParams;
       linkToAdd = linkToAdd + 1;
   end
end
fprintf('Total mass of trained model : %f \n',total_estimated_mass)

baseError = identifiableSubspaceBasis'*(inertialParametersFromModel(myModel.modelParams)-finalInertialParams);
notIdentifiableError = notIdentifiableSubspaceBasis'*(inertialParametersFromModel(myModel.modelParams)-finalInertialParams);

fprintf('Norm of the final error in base parameters : %f\n', norm(baseError))
fprintf('Norm of the final error in not identifiable space : %f\n', norm(notIdentifiableError))
identifiableError = norm(baseError);
nonIdentifiableError = norm(notIdentifiableError);


% Depending on the number of samples 'n' the updates can
% be quite relevant. With increasing 'n' the magnitude of
% the updates tends to become smaller.

for i = length(dir_ind(1:NB*6))+1 : length(dir_ind)
   cov_ini = get_field(realBNEAIP.bnt.bnet.CPD{dir_ind(i)}, 'cov');
   cov_est = get_field(bnetHat.CPD{dir_ind(i)},         'cov');
   cov_upd = cov_est - cov_ini;
   fprintf('[INFO] %s was updated by %f \n', realBNEAIP.bnt.nodes.labels{i}, norm(cov_upd)./norm(cov_ini));
end
