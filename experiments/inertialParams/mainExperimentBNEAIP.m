
res = 0;

S_dmodel  = 1e-2;
S_ymodel  = 1e-4;

NB_BNEAIP = 4;
dmodel_BNEAIP = autoTree(NB_BNEAIP);
dmodel_BNEAIP = autoTreeStochastic(dmodel_BNEAIP, S_dmodel);
dmodel_BNEAIP.gravity = [0; -9.81; 0];
nrOfSamples = 30;
nrOfTrials = 10;

% for each trial the results are stored in a matrix
% number of sensors setups  \times 2 

for trail = 1:nrOfTrials

    
nrOfSensorsSetups = 3;
initialError{trail} = zeros(1,2);
results{trail} = zeros(nrOfSensorsSetups,2);
sensorSetupTested = 1;
    
% We will test the fitting properties of several distribution of sensors

fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIAL %d ~~~~~~~~~~~~~~~~~~~~~\n',trail);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

% first we try for fitting without force sensing (but with joint accelerations) 
% cleaerly we get no improvement 
fprintf('~~~~~~~~~Testing with just joint accelerations~~~~~~~~~\n');
ymodel_BNEAIP = insufficientSensing(dmodel_BNEAIP);
ymodel_BNEAIP  = autoSensStochastic(ymodel_BNEAIP, S_ymodel);

[result{trail}(sensorSetupTested,1), result{trail}(sensorSetupTested,2), initialError{trail}(1,1), initialError{trail}(1,2) ] = experimentLearnBNEAIP(dmodel_BNEAIP, ymodel_BNEAIP, trail, nrOfSamples);
sensorSetupTested = sensorSetupTested+1;

% first we try for fitting with the base force sensing, iagnemma style 
fprintf('~~~~~~~~~Testing with base force torque + joint accelerations~~~~~~~~~\n');
ymodel_BNEAIP = basicSensing(dmodel_BNEAIP);
ymodel_BNEAIP  = autoSensStochastic(ymodel_BNEAIP, S_ymodel);

[result{trail}(sensorSetupTested,1), result{trail}(sensorSetupTested,2)] = experimentLearnBNEAIP(dmodel_BNEAIP, ymodel_BNEAIP, trail, nrOfSamples);
sensorSetupTested = sensorSetupTested+1;

% then we suppose that we know the spatial acceleration of each link 
fprintf('~~~~~~~~~Testing with base force torque + joint accelerations + link accelerations~~~~~~~~~\n');
ymodel_BNEAIP = basicSensingAndAcc(dmodel_BNEAIP);
ymodel_BNEAIP  = autoSensStochastic(ymodel_BNEAIP, S_ymodel);

[result{trail}(sensorSetupTested,1), result{trail}(sensorSetupTested,2)] = experimentLearnBNEAIP(dmodel_BNEAIP, ymodel_BNEAIP, trail, nrOfSamples);
sensorSetupTested = sensorSetupTested+1;


end