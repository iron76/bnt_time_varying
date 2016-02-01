function [ data ] = organiseBERDYCompatibleSensorData( data , subjectID, trialID)
%ORGANISEBERDYCOMPATIBLESENSORDATA Reorganises the sensor data from human
%capture into a BERDY compatible form
%   loads sensor transfors and organises data to generate the y and ys
%   matrices in correct form

%% load the processed sensor data
load('./experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat','processedSensorData');

%% load the sensor link transforms
load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');

%% Extract the time stamps of interest by looking for spikes in FT 
% Since the beginning and end of the experiment featured a little hop to
% generate a spike in the FT.

t = processedSensorData(subjectID,trialID).t;
f_temp = processedSensorData(subjectID,trialID).f_fp;
[~,chosenF_ID] = max(f_temp(:,1));
[~,tminIndex] = max(f_temp(chosenF_ID,1:round(end/2)));
[~,tmaxIndex] = max(f_temp(chosenF_ID,round(end/2):end));
totPointsInConsideration = tmaxIndex - tminIndex;
tminIndex = tminIndex + round(0.1*totPointsInConsideration);
tmaxIndex = tmaxIndex + round(length(t)/2) - round(0.1*totPointsInConsideration);

data.min_time = t(tminIndex);
data.max_time = t(tmaxIndex);
data.nsamples = tmaxIndex - tminIndex;
data.time = t(tminIndex:tmaxIndex);

data.q1 = processedSensorData(subjectID,trialID).q1(tminIndex:tmaxIndex);
data.q2 = processedSensorData(subjectID,trialID).q2(tminIndex:tmaxIndex);
data.dq1 = processedSensorData(subjectID,trialID).dq1(tminIndex:tmaxIndex);
data.dq2 = processedSensorData(subjectID,trialID).dq2(tminIndex:tmaxIndex);
data.ddq1 = processedSensorData(subjectID,trialID).ddq1(tminIndex:tmaxIndex);
data.ddq2 = processedSensorData(subjectID,trialID).ddq2(tminIndex:tmaxIndex);

data.q = [data.q1 data.q2]';
data.dq = [data.dq1 data.dq2]';
data.ddq = [data.ddq1 data.ddq2]';

%% data from IMU sensor
%Extracting the transforms
a_imu = processedSensorData(subjectID,trialID).a_imu;
fp = processedSensorData(subjectID,trialID).f_fp;

numTSteps = size(a_imu(:, (tminIndex:tmaxIndex)),2);
X_2_imu = sensorLinkTransforms.X_2_imu;
y_2_imu = X_2_imu * [zeros(3,numTSteps);processedSensorData(subjectID,trialID).a_imu(:, (tminIndex:tmaxIndex))];
y_1_fp = sensorLinkTransforms.XStar_1_fpini * fp(:, (tminIndex:tmaxIndex));


% Putting the IMU Data, accelerometer only, gyroscope consists of 0.
%data.y_imu = y_2_imu;
data.ys_imu = y_2_imu;
%data.y_ft = y_1_fp;
data.ys_fts = y_1_fp;
%data.ys_fts = data.y_fts;
save('./experiments/humanFixedBase/intermediateDataFiles/berdyFormattedSensorData.mat','data');


end
