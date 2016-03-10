% loadModelFromURDF
% Script to generate a Featherstone spatial v1/drake model from the URDF.
% The generated model is stored within a mat file. The option of
% loadingFromDrake is used to reload the mat file instead of reloading from
% drake. That way, the drake loading can be performed by users who have
% drake and spatial installed alongwith bnt_time_varying


clc; close all;clear;

%% loading subject parameters
load('./experiments/humanFixedBase/data/subjectSizeParams.mat');
addpath(genpath('./human_models/'));

subjectIDList = 1;
   
fprintf('Loading the model from URDF...\n');
fprintf('-----------------------------\n\n');
  
%adding drake path
path_to_drake_distro = '.././drake-distro/';
addpath(path_to_drake_distro);
addpath(strcat(path_to_drake_distro,'drake'));
addpath(strcat(path_to_drake_distro,'build/matlab'));
addpath_drake;

for subjectID = subjectIDList
    
    str = sprintf('./human_models/threeLinkHumanLikeSubj%02d/threeLinkHuman_subject%d.urdf',subjectID,subjectID);
    R = RigidBodyManipulator(str);
    humanThreeLink_dmodel = R.featherstone;
    humanThreeLink_dmodel.NB = 2;
    humanThreeLink_dmodel.jtype = {'R','R'};
    humanThreeLink_dmodel.appearance = {'1' '1'}';
    humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
    humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 
        
    % storing model per each subject
    humanThreeLinkModelFromURDF(subjectID).dmodel = humanThreeLink_dmodel;    
end

save('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat','humanThreeLinkModelFromURDF');
   
%removing drake path
rmpath_drake;
rmpath('.././drake-distro/');
rmpath(strcat(path_to_drake_distro,'drake'));
rmpath(strcat(path_to_drake_distro,'build/matlab'));

fprintf('\nFinished Loading the model\n');
fprintf('-----------------------------\n\n');

