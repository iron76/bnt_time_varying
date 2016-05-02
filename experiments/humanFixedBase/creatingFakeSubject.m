% %% CREATION OF SUBJECT13 : subject with URDF_2 and measurements y_3
% 
% load ('./experiments/humanFixedBase/data/imuExtractedDataGen16.mat');
% imuData(13,:) = imuData(3,:);
% save('./experiments/humanFixedBase/data/imuExtractedDataGen16_fake.mat');
%  
% load('./experiments/humanFixedBase/data/VICONsaveDataGen16.mat');
% subjectData(13,:)  =  subjectData(3,:);
% subjects(1,13)  =  17;
% trials(13,:)  =  trials(3,:);
% save('./experiments/humanFixedBase/data/VICONsaveDataGen16_fake.mat');
% 
% fprintf('\nStarting computeSubjectSpecificURDFPrams computation\n');
% subjectParams(13)  =  subjectParams(2);
% save('./experiments/humanFixedBase/data/subjectSizeParams_fake.mat','subjectParams');
% 
% %% CREATION OF SUBJECT14 : subject with URDF_5kg and measurements y_3
% 
% load ('./experiments/humanFixedBase/data/imuExtractedDataGen16_fake.mat');
% imuData(14,:) = imuData(3,:);
% save('./experiments/humanFixedBase/data/imuExtractedDataGen16_fake.mat');
% 
load('./experiments/humanFixedBase/data/VICONsaveDataGen16_fake.mat');
subjectData(14,:)  =  subjectData(3,:);
subjects(1,14)  =  18;
trials(14,:)  =  trials(3,:);
save('./experiments/humanFixedBase/data/VICONsaveDataGen16_fake.mat','subjectData','subjects','trials');
% 
% load('./experiments/humanFixedBase/data/subjectSizeParams.mat');
% fprintf('\nStarting computeSubjectSpecificURDFPrams computation\n');
% subjectParams(14)  =  subjectParams(2);
% save('./experiments/humanFixedBase/data/subjectSizeParams_fake.mat','subjectParams');

%% CREATION OF SUBJECT15 : subject with URDF_5kg and measurements y_2

load ('./experiments/humanFixedBase/data/imuExtractedDataGen16_fake.mat');
imuData(15,:) = imuData(2,:);
save('./experiments/humanFixedBase/data/imuExtractedDataGen16_fake.mat','imuData');

load('./experiments/humanFixedBase/data/VICONsaveDataGen16_fake.mat');
subjectData(15,:)  =  subjectData(2,:);
subjects(1,15)  =  19;
trials(15,:)  =  trials(2,:);
save('./experiments/humanFixedBase/data/VICONsaveDataGen16_fake.mat');

load('./experiments/humanFixedBase/data/subjectSizeParams_fake.mat');
fprintf('\nStarting computeSubjectSpecificURDFPrams computation\n');
subjectParams(15)  =  subjectParams(2);
save('./experiments/humanFixedBase/data/subjectSizeParams_fake.mat','subjectParams');