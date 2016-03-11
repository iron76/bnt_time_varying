% CREATEURDFMODELFROMSUBJECTPARAMS
% Script tp create for each subject a file urdf.

load('./experiments/humanFixedBase/data/subjectSizeParams.mat');

subjectList = 1:12;

for subjectID = subjectList
    
    urdfModelTemplate = fileread('./human_models/threeLinkHumanLikeTemplate/threeLinkHumanTemplate.urdf');
    
    urdfModelTemplate = strrep(urdfModelTemplate,'FOOTMASS',num2str(subjectParams(subjectID).footMass));
    urdfModelTemplate = strrep(urdfModelTemplate,'LEGMASS',num2str(subjectParams(subjectID).legMass));
    urdfModelTemplate = strrep(urdfModelTemplate,'TORSOMASS',num2str(subjectParams(subjectID).torsoMass));

    urdfModelTemplate = strrep(urdfModelTemplate,'FOOTALFA',num2str(subjectParams(subjectID).footWidth));
    urdfModelTemplate = strrep(urdfModelTemplate,'FOOTBETA',num2str(subjectParams(subjectID).footHeight));
    urdfModelTemplate = strrep(urdfModelTemplate,'FOOTGAMMA',num2str(subjectParams(subjectID).footDepth));

    urdfModelTemplate = strrep(urdfModelTemplate,'LEGALFA',num2str(subjectParams(subjectID).legWidth));
    urdfModelTemplate = strrep(urdfModelTemplate,'LEGBETA',num2str(subjectParams(subjectID).legHeight));
    urdfModelTemplate = strrep(urdfModelTemplate,'LEGGAMMA',num2str(subjectParams(subjectID).legDepth));

    urdfModelTemplate = strrep(urdfModelTemplate,'TORSOALFA',num2str(subjectParams(subjectID).torsoWidth));
    urdfModelTemplate = strrep(urdfModelTemplate,'TORSOBETA',num2str(subjectParams(subjectID).torsoHeight));
    urdfModelTemplate = strrep(urdfModelTemplate,'TORSOGAMMA',num2str(subjectParams(subjectID).torsoDepth));

    urdfModelTemplate = strrep(urdfModelTemplate,'HALFfootBETA',num2str(0.5*subjectParams(subjectID).footHeight));
    urdfModelTemplate = strrep(urdfModelTemplate,'HALFlegBETA',num2str(0.5*subjectParams(subjectID).legHeight));
    urdfModelTemplate = strrep(urdfModelTemplate,'HALFtorsoBETA',num2str(0.5*subjectParams(subjectID).torsoHeight));

    urdfModelTemplate = strrep(urdfModelTemplate,'FOOTINERTIAIXX',num2str(subjectParams(subjectID).footIxx));
    urdfModelTemplate = strrep(urdfModelTemplate,'FOOTINERTIAIYY',num2str(subjectParams(subjectID).footIyy));
    urdfModelTemplate = strrep(urdfModelTemplate,'FOOTINERTIAIZZ',num2str(subjectParams(subjectID).footIzz));

    urdfModelTemplate = strrep(urdfModelTemplate,'LEGINERTIAIXX',num2str(subjectParams(subjectID).legIxx));
    urdfModelTemplate = strrep(urdfModelTemplate,'LEGINERTIAIYY',num2str(subjectParams(subjectID).legIyy));
    urdfModelTemplate = strrep(urdfModelTemplate,'LEGINERTIAIZZ',num2str(subjectParams(subjectID).legIzz));

    urdfModelTemplate = strrep(urdfModelTemplate,'TORSOINERTIAIXX',num2str(subjectParams(subjectID).torsoIxx));
    urdfModelTemplate = strrep(urdfModelTemplate,'TORSOINERTIAIYY',num2str(subjectParams(subjectID).torsoIyy));
    urdfModelTemplate = strrep(urdfModelTemplate,'TORSOINERTIAIZZ',num2str(subjectParams(subjectID).torsoIzz));

    urdfModelTemplate = strrep(urdfModelTemplate,'ANKLEVAL',num2str((-subjectParams(subjectID).footDepth)/4));

    %% storing file into a proper folder

    fileID = fopen(sprintf('threeLinkHuman_subject%d.urdf',subjectID),'w');
    fprintf(fileID,'%s', urdfModelTemplate);
    fclose(fileID);

    movefile(sprintf('threeLinkHuman_subject%d.urdf',subjectID),sprintf('./human_models/threeLinkHumanLikeSubj%02d',subjectID));

end
fprintf('---------\n');
fprintf('URDF models successfully created!\n');
