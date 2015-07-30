
res = 0;

S_dmodel  = 1e-2;
S_ymodel  = 1e-4;

NB_BNEAIP = 4;
dmodel_BNEAIP = autoTree(NB_BNEAIP);
dmodel_BNEAIP = autoTreeStochastic(dmodel_BNEAIP, S_dmodel);
dmodel_BNEAIP.gravity = [0; -9.81; 0];
nrOfSamples = 30;
nrOfTrials = 10;

for trail = 1:nrOfTrials

% We will test the fitting properties of several distribution of sensors

fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIAL %d ~~~~~~~~~~~~~~~~~~~~~\n',trail);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

% first we try for fitting without force sensing (but with joint accelerations) 
% cleaerly we get no improvement 
fprintf('~~~~~~~~~Testing with just joint accelerations~~~~~~~~~\n');
ymodel_BNEAIP = insufficientSensing(dmodel_BNEAIP);
ymodel_BNEAIP  = autoSensStochastic(ymodel_BNEAIP, S_ymodel);

experimentLearnBNEAIP(dmodel_BNEAIP, ymodel_BNEAIP, trail, nrOfSamples);

% first we try for fitting with the base force sensing, iagnemma style 
fprintf('~~~~~~~~~Testing with base force torque + joint accelerations~~~~~~~~~\n');
ymodel_BNEAIP = basicSensing(dmodel_BNEAIP);
ymodel_BNEAIP  = autoSensStochastic(ymodel_BNEAIP, S_ymodel);

experimentLearnBNEAIP(dmodel_BNEAIP, ymodel_BNEAIP, trail, nrOfSamples);

% then we suppose that we know the spatial acceleration of each link 
fprintf('~~~~~~~~~Testing with base force torque + joint accelerations + link accelerations~~~~~~~~~\n');
ymodel_BNEAIP = basicSensingAndAcc(dmodel_BNEAIP);
ymodel_BNEAIP  = autoSensStochastic(ymodel_BNEAIP, S_ymodel);

experimentLearnBNEAIP(dmodel_BNEAIP, ymodel_BNEAIP, trail, nrOfSamples);

end