
path_to_drake_distro = '/home/naveenoid/Workspace/Simulation/Drake/drake-distro/';
addpath(path_to_drake_distro);
addpath(strcat(path_to_drake_distro,'drake'));
addpath(strcat(path_to_drake_distro,'build/matlab'));

addpath_drake;

%
R = RigidBodyManipulator('./robot_models/fiveLinkHumanLike/fiveLinkHuman.urdf');
humanFiveLink_dmodel = R.featherstone;
humanFiveLink_dmodel.NB = 4;

humanFiveLink_dmodel.jtype = {'R','R','R','R'};
humanFiveLink_dmodel.appearance = {'1' '1','1','1'}';

rmpath_drake;

rmpath(path_to_drake_distro);
rmpath(strcat(path_to_drake_distro,'drake'));
rmpath(strcat(path_to_drake_distro,'build/matlab'));
