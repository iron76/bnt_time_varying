function [ data ] = organiseBERDYCompatibleSensorData( data , subjectID, trialID)
%ORGANISEBERDYCOMPATIBLESENSORDATA Reorganises the sensor data from human
%capture into a BERDY compatible form
%   loads sensor transfors and organises data to generate the y and ys
%   matrices in correct form


load('./experiments/humanFixedBase/data/preProcessedSensorData.mat','processedSensorData');

t = processedSensorData(subjectID,trialID).t;
f_temp = processedSensorData(subjectID,trialID).f_0;
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

%% data from force plate sensing
%f_0 from sensor frame extraction is angular-linear. The MAP y requires
%linear-angular fts measurement

data.y_fts = zeros(size(processedSensorData(subjectID,trialID).f_0(:,(tminIndex:tmaxIndex))));
data.y_fts(1:3,:) = - processedSensorData(subjectID,trialID).f_0(4:6, (tminIndex:tmaxIndex));
data.y_fts(4:6,:) = - processedSensorData(subjectID,trialID).f_0(1:3, (tminIndex:tmaxIndex));

data.ys_fts = data.y_fts;


%% data from IMU sensor

%IMU from sensor frame extraction is angular-linear. The MAP y requires
%linear-angular IMU measurement
data.y_imu = [processedSensorData(subjectID,trialID).a_2_imulin(:, (tminIndex:tmaxIndex));
              zeros(3,length(tminIndex:tmaxIndex))];
% data.y_imu = [zeros(3,length(tminIndex:tmaxIndex));
%               processedSensorData(subjectID,trialID).a_2_imulin(:, (tminIndex:tmaxIndex))];
data.ys_imu = data.y_imu;


save('./experiments/humanFixedBase/data/processedSensorData.mat','data');


end
