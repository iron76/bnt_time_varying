if(exist('./experiments/humanFixedBase/humanThreeLinkModelFromURDF.mat','file'))
    fprintf('Loading preconfigured model\n');
    fprintf('-----------------------------\n\n');
    load('./experiments/humanFixedBase/humanThreeLinkModelFromURDF.mat');
else

    fprintf('Loading the model from URDF\n');
    fprintf('-----------------------------\n\n');
%     
%     path_to_drake_distro = '/drake-distro/';
%     addpath(path_to_drake_distro);
%     addpath(strcat(path_to_drake_distro,'drake'));
%     addpath(strcat(path_to_drake_distro,'build/matlab'));
% 
%     addpath_drake;

    %
    str = sprintf('./robot_models/threeLinkHumanLikeSubj%d/threeLinkHuman_subject%d.urdf',subjectID,subjectID)
    R = RigidBodyManipulator(sprintf('./robot_models/threeLinkHumanLikeSubj%d/threeLinkHuman_subject%d.urdf',subjectID,subjectID));
    humanThreeLink_dmodel = R.featherstone;
    humanThreeLink_dmodel.NB = 2;

    humanThreeLink_dmodel.jtype = {'R','R'};
    humanThreeLink_dmodel.appearance = {'1' '1'}';
    humanThreeLink_dmodel.linkname = {'foot' 'leg' 'torso'}; 

%     rmpath_drake;
% 
%     rmpath(path_to_drake_distro);
%     rmpath(strcat(path_to_drake_distro,'drake'));
%     rmpath(strcat(path_to_drake_distro,'build/matlab'));
%     humanThreeLink

    
    save('./experiments/humanFixedBase/humanThreeLinkModelFromURDF.mat','humanThreeLink_dmodel');
end
fprintf('\nFinished Loading the model\n');
fprintf('-----------------------------\n\n');