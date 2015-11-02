% if(exist(sprintf('./experiments/humanFixedBase/humanThreeLinkModelFromURDF_subject%d.mat',subjectID),'file'))
%     fprintf('Loading preconfigured model...\n');
%     %fprintf('-----------------------------\n\n');
%     load(sprintf('./experiments/humanFixedBase/humanThreeLinkModelFromURDF_subject%d.mat',subjectID));
% else
%     
%     

loadingFromDrake = 0;
if (loadingFromDrake)
    fprintf('Loading the model from URDF...\n');
    %fprintf('-----------------------------\n\n');
    
    path_to_drake_distro = '../../Drake/drake-distro/';
    addpath(path_to_drake_distro);
    addpath(strcat(path_to_drake_distro,'drake'));
    addpath(strcat(path_to_drake_distro,'build/matlab'));

    addpath_drake;

    for subjID = 1:3
   
        str = sprintf('./robot_models/threeLinkHumanLikeSubj%d/threeLinkHuman_subject%d.urdf',subjID,subjID)
        R = RigidBodyManipulator(str);
        humanThreeLink_dmodel = R.featherstone;
        humanThreeLink_dmodel.NB = 2;

        humanThreeLink_dmodel.jtype = {'R','R'};
        humanThreeLink_dmodel.appearance = {'1' '1'}';
        humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
        humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 
        
        save(sprintf('./experiments/humanFixedBase/humanThreeLinkModelFromURDF_subject%d.mat',subjID),'humanThreeLink_dmodel');
    end
     rmpath_drake;
   
  
    rmpath(path_to_drake_distro);
    rmpath(strcat(path_to_drake_distro,'drake'));
    rmpath(strcat(path_to_drake_distro,'build/matlab'));

end
 load(sprintf('./experiments/humanFixedBase/humanThreeLinkModelFromURDF_subject%d.mat',subjID));

fprintf('\nFinished Loading the model\n');
fprintf('-----------------------------\n\n');