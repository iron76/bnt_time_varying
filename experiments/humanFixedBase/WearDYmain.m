
clear;clc;close all;
% forcing the casual generator to be const
rng(1);

%% testOptions
plotTestQ = true;
plotLogLandQ = true;
plotSigma = false;

%% add\create folders

dataFolder = './experiments/humanFixedBase/data/';
if (exist (dataFolder,'dir')==7)
    addpath(genpath(dataFolder));
end

helperFunctionsFolder = './experiments/humanFixedBase/helperFunctions/';
if (exist (helperFunctionsFolder,'dir')==7)
    addpath(genpath(helperFunctionsFolder));
end

intermediateDataFilesFolder = './experiments/humanFixedBase/intermediateDataFiles/';
if (exist(intermediateDataFilesFolder,'dir')==7)
    addpath(genpath(intermediateDataFilesFolder));
end

testsFolder = './experiments/humanFixedBase/tests/';
if (exist(testsFolder,'dir')==7)
    addpath(genpath(testsFolder));
end

%% add\create files .mat

fprintf('Starting synchronisedData computation\n');
file = './experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat';
if(exist(file,'file')==2)
    load(file);
else    
    synchroniseCaptureData
end
%
fprintf('\nStarting computeSubjectSpecificURDFPrams computation\n');
file = './experiments/humanFixedBase/data/subjectSizeParams.mat';
if(exist(file,'file')==2)
    load(file);
else    
    computeSubjectSpecificURDFParams
end
%
fprintf('\nStarting createUrdfModelFromSubjectParam computation\n');
folder = './human_models/';
if(exist(folder,'dir')==7)
    addpath(genpath(folder))
else    
    createUrdfModelFromSubjectParams
end
%
fprintf('\nStarting loadModelFromURDF computation\n');  
    loadModelFromURDF
%
fprintf('\nStarting computeLinkSensorFrames computation\n');
file = './experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat';
if(exist(file,'file')==2)
    load(file);
else    
    computeLinkSensorFrames
end
%
fprintf('\nStarting organiseBERDYCompatibleSensorData computation\n');
file = './experiments/humanFixedBase/intermediateDataFiles/BERDYFormattedSensorData.mat';
if(exist(file,'file')==2)
    load(file);
else    
    organiseBERDYCompatibleSensorData
end

%% selected subjects and trials
subjectList = 1;
trialList = 1;  

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d ',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d\n',trialID);
        
        %% load BERDY compatible data and sensor link transforms
        load('./experiments/humanFixedBase/intermediateDataFiles/BERDYFormattedSensorData.mat');

        currentTrial = BERDYFormattedSensorData(subjectID,trialID); 
        data = currentTrial.data;
        dataTime = data.dataTime;
        q = currentTrial.data.q';
        dq = currentTrial.data.dq';
        ddq = currentTrial.data.ddq';
        len = length(dataTime);

        load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');

        currentTrialSens = sensorLinkTransforms(subjectID,trialID);
        X_imu_2 = currentTrialSens.X_imu_2;
        XStar_fp_0 = currentTrialSens.XStar_fp_0;
        XStar_0_1 = currentTrialSens.XStar_0_1;

        %% build model valid for all MAP case 
        load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');

        currentModel = humanThreeLinkModelFromURDF(subjectID);

        humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
        humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 
        
        dmodel  = currentModel.dmodel;                      % deterministic model
        dmodel  = autoTreeStochastic(dmodel, 1e1);         % probabilistic model for D equation (added Sv and Sw)
        myModel = model(dmodel);    

        load('./experiments/humanFixedBase/data/subjectSizeParams.mat');
        currentParams = subjectParams(subjectID);
        
        
        [myMAP, Ymatrix,b_Y,data.y] = setSensorsNumber(dmodel, myModel, len,q,dq,ddq,currentParams, currentTrialSens,data);
        
        %% Test 
        %  to verify the concavity of Q function wrt a diagonal
        %  element variation of theta 
% %     
% %         len = 1;
% %         SigmaTest = eye(20);
% % 
% %         params= linspace(1e-6, 1e2,100);
% % %         params_d2q2 = linspace(1e-3, 0.5,100);
% %         l_test = zeros(length(params),2);
% %                 
% %         for i = 1:length(params)
% %             SigmaTest(2,2) = params(i);
% %             for j = 1: length(params_d2q)
% %                 
% %                 SigmaTest(1,1) = params(j);
% %                 l_test(j,i) =0;
% %                 for t = 1:len
% %                 [l_tmp] = computeLogL(myMAP,q(:,t),dq(:,t),Ymatrix{t}, data.y(:,t),b_Y(:,t), SigmaTest);
% %                 l_test(j,i) = l_test(j,i) + l_tmp;
% %                 end
% %             end
% %         end
%
% % % %         params_d2q1 = linspace(1e-3, 0.5,100);
% % % %         params_d2q2 = linspace(1e-3, 0.5,100);
% % % %         l_test = zeros(length(params_d2q1)*length(params_d2q1),3);
% % % %         
% % % %         
% % % %         z = 1;
% % % %         for i = 1:length(params_d2q2)
% % % %             SigmaTest(2,2) = params_d2q2(i);
% % % %             for j = 1: length(params_d2q1)
% % % %                 len = 1;
% % % %                 SigmaTest(1,1) = params_d2q1(j);
% % % %                 [l_tmp] = computeLogL(myMAP_1sens,q(:,len),dq(:,len),Ymatrix_1sens{len}, data.y_1sens(:,len),b_Y_1sens(:,len), SigmaTest);
% % % %                 l_test(z, :) = [SigmaTest(1,1), SigmaTest(2,2),  l_tmp];
% % % %                 
% % % %                 z = z + 1;
% % % %             end
% % % %         end

          %% Computing EM
         
        fprintf('\nStarting EM procedure\n');    
        
%         vect = [40;50];
%         SigmaTest = diag(vect);
        
%         sigma_ygivend_init = SigmaTest;      
        sigma_ygivend_init = inv(full(myMAP.IDsens.sensorsParams.Sy_inv));

        
%         [Sigma_EM,l_EM] = computeEM_test(myMAP_4sens,data.q',data.dq' ,Ymatrix_4sens, data.y_4sens,b_Y_4sens, sigma_ygivend_init,'MAX_ITER',200);
          [Sigma_EM, Qin_EM, lin_EM, Qfin_EM, lfin_EM] = computeEM_test(myMAP,data.q',data.dq' ,Ymatrix, data.y,b_Y, sigma_ygivend_init,'MAX_ITER',200);
%          [Sigma_EM, l_EM] = computeEM_test(myMAP_1sens,data.q',data.dq' ,Ymatrix_1sens, data.y_1sens,b_Y_1sens, sigma_ygivend_init,'MAX_ITER',200);
 
        
         
         
         
         
        %% Inequalities plot
        
        fig = figure();
        axes1 = axes('Parent',fig,'FontSize',16);
        box(axes1,'on');
        hold(axes1,'on');

        plot1 = plot(lin_EM(:,1), 'lineWidth',2.5,'LineStyle','-'); hold on;
        plot2 = plot(Qin_EM, 'lineWidth',2.5,'LineStyle','-'); hold on;
        plot3 = plot(Qfin_EM, 'lineWidth',2.5,'LineStyle','-'); hold on;
        plot4 = plot(lfin_EM(:,1), 'lineWidth',2.5,'LineStyle','-'); hold on;
        
        leg = legend('$\ell(\theta^k)$','$\mathcal Q(\theta^k,\theta^k)$', '$\mathcal Q(\theta^k,\theta^{k+1})$','$\ell(\theta^{k+1})$','Location','northeast');
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel(' k iterations','FontSize',20);
        ylabel('inequalities','FontSize',20);
        axis tight;
        

        
% % %         %% Plot log-likelihood
% % %         
% % % % %         % folder for plots
% % % % %         figFolder = './experiments/humanFixedBase/plots';
% % % % %         if(exist(figFolder,'dir')==0)
% % % % %             mkdir(figFolder);
% % % % %         end
% % % % %         
% % % % %         fig = figure();
% % % % %         axes1 = axes('Parent',fig,'FontSize',16);
% % % % %         box(axes1,'on');
% % % % %         hold(axes1,'on');
% % % % %         %grid on;
% % % % %         
% % % % %         plot1 = plot(l_EM, 'lineWidth',2.5,'LineStyle','-'); hold on;
% % % % %         set(plot1,'color',[0 0 0]);
% % % % %         %xlabel(' k iterations','FontSize',20);
% % % % %         ylabel('log-likelihood','FontSize',20);
% % % % %         axis tight;
% % % % %         title('Parameter $\theta = \Sigma_{y|d}^{-1}$',...
% % % % %         'HorizontalAlignment','center',...
% % % % %         'FontSize',20,...
% % % % %         'Interpreter','latex');
% % % % %         %grid on;
% % % % %         
% % %         
% % % %%
% % %         %% Plot log-likelihood and auxiliary function Q
% % %         
% % %         % folder for plots
% % % % %         figFolder = './experiments/humanFixedBase/plots';
% % % % %         if(exist(figFolder,'dir')==0)
% % % % %             mkdir(figFolder);
% % % % %         end
% % % % %         
% % % % %         fig = figure();
% % % % %         axes1 = axes('Parent',fig,'FontSize',16);
% % % % %         box(axes1,'on');
% % % % %         hold(axes1,'on');
% % % % %         %grid on;
% % % % %         
% % % % %         plot1 = plot(Q, 'lineWidth',2.5,'LineStyle','-'); hold on;
% % % % %         plot2 = plot(logL, 'lineWidth',2.5,'LineStyle','-'); 
% % % % %         %set(plot1,'color',[0 0 0]);
% % % % %         leg = legend('$\mathcal Q$','$\ell$','Location','southwest');
% % % % %                 set(leg,'Interpreter','latex');
% % % % %                 set(leg,'FontSize',18);
% % % % %         xlabel(' k iterations','FontSize',20);
% % % % %         %ylabel('Q','FontSize',20);
% % % % %         axis tight;
% % % % % 
% % % % %         
% % % % %         fig = figure();
% % % % %         subplot (211)
% % % % %         plot1 = plot(logL, 'lineWidth',2.5,'LineStyle','-'); hold on;
% % % % %         set(plot1,'color',[0 0 0]);
% % % % %         %xlabel(' k iterations','FontSize',20);
% % % % %         ylabel('$\ell$','FontSize',20,...
% % % % %         'Interpreter','latex');
% % % % %         axis tight;
% % % % %         title('Parameter $\theta = \Sigma_{y|d}^{-1}$',...
% % % % %         'HorizontalAlignment','center',...
% % % % %         'FontSize',20,...
% % % % %         'Interpreter','latex');
% % % % %         %grid on;
% % % % %         
% % % % %         subplot (212)
% % % % %         plot1 = plot(Q, 'lineWidth',2.5,'LineStyle','-'); hold on;
% % % % %         set(plot1,'color',[0 0 0]);
% % % % %         %xlabel(' k iterations','FontSize',20);
% % % % %         ylabel('$\mathcal Q$','FontSize',20,...
% % % % %         'Interpreter','latex');
% % % 
% % %         %% Plot concavity of Q
% % % % % %         if (plotTestQ)
% % % % % %             fig = figure();
% % % % % %             axes1 = axes('Parent',fig,'FontSize',16);
% % % % % %             box(axes1,'on');
% % % % % %             hold(axes1,'on');
% % % % % %             %grid on;
% % % % % % 
% % % % % %             %legendInfo = cell{1, length(Q_paramTest)};
% % % % % %             for lenIter = 1:length(Q_paramTest)
% % % % % %                 plot1 = plot(Q_paramTest{lenIter}, 'lineWidth',2.5,'LineStyle','-'); hold on;
% % % % % % 
% % % % % %                 legendInfo{lenIter}=['$\mathcal Q$', num2str(lenIter)]; 
% % % % % %                 leg = legend(legendInfo)
% % % % % %                     set(leg,'Interpreter','latex');
% % % % % %                     set(leg,'FontSize',15);
% % % % % %             end
% % % % % %         end
% % % % %         %%
% % % % %         if(plotLogLandQ)
% % % % %         save2pdf(fullfile(figFolder, (sprintf('Subject %d,Trial %d, logL',subjectID,trialID))),fig,600);
% % % % %         end
% % % % %         
% % % % %         if(plotSigma)
% % % % %         for k = 1:length(logL)
% % % % %            clims = [-5,5];  %limits of the colorbar
% % % % %            fig = figure();
% % % % %            imagesc(Sigma_EM{k, 1},clims)
% % % % %            colorbar
% % % % %            title(sprintf('Covariance at iteration %d',k))
% % % % %         end
% % % % %         end

    end
    fprintf('\n');
end
%% storing results
%save('./experiments/humanFixedBase/intermediateDataFiles/finalResults.mat','finalResults');

fprintf('---------\n');
fprintf('Done!\n');

