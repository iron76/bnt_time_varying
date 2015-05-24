fprintf('Loading the model from URDF\n');
fprintf('-----------------------------\n\n');
path_to_drake_distro = '/home/naveenoid/Workspace/Simulation/Drake/drake-distro/';
addpath(path_to_drake_distro);
addpath(strcat(path_to_drake_distro,'drake'));
addpath(strcat(path_to_drake_distro,'build/matlab'));

addpath_drake;

%
R = RigidBodyManipulator(sprintf('./robot_models/threeLinkHumanLike/threeLinkHuman_subject%d.urdf',subjectID));
humanThreeLink_dmodel = R.featherstone;
humanThreeLink_dmodel.NB = 2;

humanThreeLink_dmodel.jtype = {'R','R'};
humanThreeLink_dmodel.appearance = {'1' '1'}';
humanThreeLink_dmodel.linkname = {'foot' 'leg' 'torso'}; 

rmpath_drake;

rmpath(path_to_drake_distro);
rmpath(strcat(path_to_drake_distro,'drake'));
rmpath(strcat(path_to_drake_distro,'build/matlab'));

fprintf('\nFinished Loading the model\n');
fprintf('-----------------------------\n\n');