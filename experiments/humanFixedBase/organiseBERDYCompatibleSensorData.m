function [ data ] = organiseBERDYCompatibleSensorData( data )
%ORGANISEBERDYCOMPATIBLESENSORDATA Reorganises the sensor data from human
%capture into a BERDY compatible form
%   loads sensor transfors and organises data to generate the y and ys
%   matrices in correct form


% if (isempty(data))
%     data=struct;
% end

load('./experiments/humanFixedBase/preProcessedSensorData.mat','processedSensorData');

subjectID = 1;
trialID = 1;
% computing min and max time by looking at peaks on f/t
t = processedSensorData(subjectID,trialID).t;
f_temp = processedSensorData(subjectID,trialID).ftx;
[val,chosenF_ID] = max(f_temp(:,1));
[val,tminIndex] = max(f_temp(chosenF_ID,1:round(end/2)));
[val,tmaxIndex] = max(f_temp(chosenF_ID,round(end/2):end));
totPointsInConsideration = tmaxIndex - tminIndex;
tminIndex = tminIndex + round(0.1*totPointsInConsideration);
tmaxIndex = tmaxIndex + round(length(t)/2) - round(0.1*totPointsInConsideration);
figure;plot(t,f_temp); axis tight;xlabel('time (sec');ylabel('Force (N)');
fprintf('tmin : %3.3f, tmax = %3.3f\n',t(tminIndex),t(tmaxIndex));

data.min_time = t(tminIndex);
data.max_time = t(tmaxIndex);
data.nsamples = tmaxIndex - tminIndex;
data.time = t(tminIndex:tmaxIndex);
% 
% data.q1 = processedSensorData(subjectID,trialID).q1;
% data.q2 = processedSensorData(subjectID,trialID).q2;
% data.dq1 = processedSensorData(subjectID,trialID).dq1;
% data.dq2 = processedSensorData(subjectID,trialID).dq2;
% data.ddq1 = processedSensorData(subjectID,trialID).ddq1;
% data.ddq2 = processedSensorData(subjectID,trialID).ddq2;
% 
% data.d2q = [data.ddq1 data.ddq2]';
% data.q_l = data.q1; data.dq_l = data.dq1; data.ddq_l = data.ddq1;
% data.q_t = data.q2; data.dq_t = data.dq2; data.ddq_t = data.ddq2;
% 
% data.y_ftx = processedSensorData(subjectID,trialID).ftx;
% data.y_fts = data.y_ftx;
% data.y_imu = processedSensorData(subjectID,trialID).imu;
% 
% data.ys_fts = data.y_fts;
% data.ys_fts = data.y_ftx;
% data.ys_ftx = data.y_ftx;
% data.ys_imu = data.y_imu;

data.q1 = processedSensorData(subjectID,trialID).q1(tminIndex:tmaxIndex);
data.q2 = processedSensorData(subjectID,trialID).q2(tminIndex:tmaxIndex);
data.dq1 = processedSensorData(subjectID,trialID).dq1(tminIndex:tmaxIndex);
data.dq2 = processedSensorData(subjectID,trialID).dq2(tminIndex:tmaxIndex);
data.ddq1 = processedSensorData(subjectID,trialID).ddq1(tminIndex:tmaxIndex);
data.ddq2 = processedSensorData(subjectID,trialID).ddq2(tminIndex:tmaxIndex);

data.d2q = [data.ddq1 data.ddq2]';
data.q_l = data.q1; data.dq_l = data.dq1; data.ddq_l = data.ddq1;
data.q_t = data.q2; data.dq_t = data.dq2; data.ddq_t = data.ddq2;

data.q = [data.q1 data.q2]';
data.dq = [data.dq1 data.dq2]';

data.y_ftx = processedSensorData(subjectID,trialID).ftx(:,(tminIndex:tmaxIndex));
data.y_fts = data.y_ftx;
data.y_imu = processedSensorData(subjectID,trialID).imu(:,(tminIndex:tmaxIndex));

data.ys_fts = data.y_fts;
data.ys_fts = data.y_ftx;
data.ys_ftx = data.y_ftx;
data.ys_imu = data.y_imu;

save('./experiments/humanFixedBase/processedSensorData.mat','data');


end
