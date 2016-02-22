function [ output_args ] = createUrdfModelFromSubjectParams( templateName)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


load('./experiments/humanFixedBase/data/subjectSizeParams.mat');

urdfModelTemplate = fileread('threeLinkHumanTemplate.urdf');

for i= 1
    
urdfModelTemplate = strrep(urdfModelTemplate,'FOOTMASS',num2str(subjectParams(i).footMass));
urdfModelTemplate = strrep(urdfModelTemplate,'LEGMASS',num2str(subjectParams(i).legMass));
urdfModelTemplate = strrep(urdfModelTemplate,'TORSOMASS',num2str(subjectParams(i).torsoMass));

urdfModelTemplate = strrep(urdfModelTemplate,'FOOTALFA',num2str(subjectParams(i).footWidth));
urdfModelTemplate = strrep(urdfModelTemplate,'FOOTBETA',num2str(subjectParams(i).footHeight));
urdfModelTemplate = strrep(urdfModelTemplate,'FOOTGAMMA',num2str(subjectParams(i).footDepth));

urdfModelTemplate = strrep(urdfModelTemplate,'LEGALFA',num2str(subjectParams(i).legWidth));
urdfModelTemplate = strrep(urdfModelTemplate,'LEGBETA',num2str(subjectParams(i).legHeight));
urdfModelTemplate = strrep(urdfModelTemplate,'LEGGAMMA',num2str(subjectParams(i).legDepth));

urdfModelTemplate = strrep(urdfModelTemplate,'TORSOALFA',num2str(subjectParams(i).torsoWidth));
urdfModelTemplate = strrep(urdfModelTemplate,'TORSOBETA',num2str(subjectParams(i).torsoHeight));
urdfModelTemplate = strrep(urdfModelTemplate,'TORSOGAMMA',num2str(subjectParams(i).torsoDepth));

urdfModelTemplate = strrep(urdfModelTemplate,'HALFfootBETA',num2str(0.5*subjectParams(i).footHeight));
urdfModelTemplate = strrep(urdfModelTemplate,'HALFlegBETA',num2str(0.5*subjectParams(i).legHeight));
urdfModelTemplate = strrep(urdfModelTemplate,'HALFtorsoBETA',num2str(0.5*subjectParams(i).torsoHeight));

urdfModelTemplate = strrep(urdfModelTemplate,'FOOTINERTIAIXX',num2str(subjectParams(i).footIxx));
urdfModelTemplate = strrep(urdfModelTemplate,'FOOTINERTIAIYY',num2str(subjectParams(i).footIyy));
urdfModelTemplate = strrep(urdfModelTemplate,'FOOTINERTIAIZZ',num2str(subjectParams(i).footIzz));

urdfModelTemplate = strrep(urdfModelTemplate,'LEGINERTIAIXX',num2str(subjectParams(i).legIxx));
urdfModelTemplate = strrep(urdfModelTemplate,'LEGINERTIAIYY',num2str(subjectParams(i).legIyy));
urdfModelTemplate = strrep(urdfModelTemplate,'LEGINERTIAIZZ',num2str(subjectParams(i).legIzz));

urdfModelTemplate = strrep(urdfModelTemplate,'TORSOINERTIAIXX',num2str(subjectParams(i).torsoIxx));
urdfModelTemplate = strrep(urdfModelTemplate,'TORSOINERTIAIYY',num2str(subjectParams(i).torsoIyy));
urdfModelTemplate = strrep(urdfModelTemplate,'TORSOINERTIAIZZ',num2str(subjectParams(i).torsoIzz));

urdfModelTemplate = strrep(urdfModelTemplate,'ANKLEVAL',num2str((-subjectParams(i).footDepth)/4));

%create subject directories
%mkdir('human_models',sprintf('threeLinkHumanLikeSubj%02d',i))


% 'human_models/threeLinkHumanLikeSubj%02d/threeLinkHuman_subject%02d','w',i,i

% %storing file into a proper folder
fileID = fopen(sprintf('threeLinkHuman_subject%d.urdf',i),'w');
fprintf(fileID,'%s', urdfModelTemplate);
fclose(fileID);

end
end

