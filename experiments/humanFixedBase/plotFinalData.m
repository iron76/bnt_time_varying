clear;clc;close all;

%% testOptions
Section1 = true;         % Section 1 plots for each subject the mean of trials/case all sensors.
    plotPaper1 = true;

Section2 = true;         % Section 2 plots q,dq,tau_MAP(all sensors).
    plotPaper2 = true;

Section3 = true;         % Section 3 plots tau and covariances adding progressively sensors
plotMAP2torques = true;
    plotPaperTau1 = true;
    plotPaperTau2 = true;
    
plotCovariances = true;
    plotPaperCov = true;
plotBarCovariances = true;
    plotPaperBarCov =true;
    
Section4 = false;    
    plotSensor =true;

%% folder for plots
figFolder = './experiments/humanFixedBase/plots';
if(exist(figFolder,'dir')==0)
    mkdir(figFolder);
end
%% ===============================SECTION1===============================%%
if(Section1)
        %% load data
        load('./experiments/humanFixedBase/intermediateDataFiles/finalResults.mat');
        

        %% cut each trial up to lenMin sample
        subjectNum =1;
        trialNum = 2;
        
        
        len = zeros(trialNum,subjectNum);

        for i = 1:subjectNum
            for j = 1:trialNum
                len(j,i) = length(finalResults(i,j).resMAP_4sens.tau_hip);
            end
        end

        lenMin = min(len);
        lenMinIndex = zeros(size(lenMin));
        for i = 1: subjectNum
            lenMinIndex(:,i) = find(len(:,i)==lenMin(:,i)); % looking for the index of the smaller sample
        end

        %% plots
        subjectList = 1;
        trialList = 1:2;

        for subjectID = subjectList
            fprintf('\n---------\nSubject : %d\n',subjectID);

            currentLenMin = lenMin(:,subjectID);

            tau_hip = zeros (trialNum,currentLenMin);
            tau_hipMean = zeros(1,currentLenMin);
            tau_ankle = zeros (trialNum,currentLenMin);
            tau_ankleMean = zeros(1,currentLenMin);

            for trialID = trialList
                fprintf('Trial : %d\n',trialID);

                currentSynchTrial = finalResults(subjectID,trialID);

                % hip
                tau_hip(trialID,:) = (currentSynchTrial.resMAP_4sens.tau_hip(1,1:currentLenMin)); 
                % ankle
                tau_ankle(trialID,:) = (currentSynchTrial.resMAP_4sens.tau_ankle(1,1:currentLenMin));
            end

            fprintf('\n');
            %% mean

            % hip 
            for i = 1:currentLenMin
                tau_hipMean(:,i) = mean(tau_hip(:,i));
            end
            
            % ankle 
            for i = 1:currentLenMin
                tau_ankleMean(:,i) = mean(tau_ankle(:,i));
            end


            %% plot for each subject comparing taus

                dataTime = finalResults(subjectID,lenMinIndex(:,subjectID)).dataTime;

                fig = figure();
                axes1 = axes('Parent',fig,'FontSize',16);
                box(axes1,'on');
                hold(axes1,'on');
                grid on;
                title(sprintf('Subject: %d, mean of trials',subjectID))

                plot1 = plot(dataTime,tau_hipMean, 'lineWidth',1.5,'LineStyle','-'); hold on;
                set(plot1,'color',[1 0 0]);
                plot2 = plot(dataTime,tau_ankleMean, 'lineWidth',1.5,'LineStyle','-'); hold on;
                set(plot2,'color',[0 0.498039215803146 0]);

                leg = legend('$\tau_{1,MAP}$','$\tau_{2,MAP}$','Location','southeast');
                set(leg,'Interpreter','latex');
                set(leg,'FontSize',18);
                xlabel('Time [s]','FontSize',20);
                ylabel('Joint torque [Nm]','FontSize',20);
                axis tight;
                grid on;

                if(plotPaper1)
                   save2pdf(fullfile(figFolder, (sprintf('Subject %d, mean of trials',subjectID))),fig,600);
                end
        end
end


%% ===============================SECTION2===============================%%
if(Section2)
        %% load data
        load('./experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat');   
        load('./experiments/humanFixedBase/intermediateDataFiles/finalResults.mat');

        subjectList = 1;
        trialList = 2;
        
        for subjectID = subjectList
            fprintf('\n---------\nSubject : %d ',subjectID);
            for trialID = trialList
                fprintf('\nTrial : %d ',trialID);

                currentSynchTrial = synchronisedData(subjectID,trialID);
                dataTime = currentSynchTrial.dataTime;
                q1 =  currentSynchTrial.q(:,1);
                q2 =  currentSynchTrial.q(:,2);
                dq1 =  currentSynchTrial.dq(:,1);
                dq2 =  currentSynchTrial.dq(:,2);
                
                currentResTrial = finalResults(subjectID,trialID);
                tau_ankle = currentResTrial.resMAP_4sens.tau_ankle;
                tau_hip = currentResTrial.resMAP_4sens.tau_hip;

                %% plots
                fig = figure();
                axes1 = axes('Parent',fig,'FontSize',16);
                box(axes1,'on');
                hold(axes1,'on');
                grid on;

                subplot(311);
                plot1 = plot(dataTime,q1.*(180/pi),'lineWidth',1.5); hold on;
                set(plot1,'color',[1 0 0]);
                plot2= plot(dataTime,q2.*(180/pi),'lineWidth',1.5); hold on;
                set(plot2,'color',[0 0.498039215803146 0]);
                leg = legend('$q_1$','$q_2$','Location','northeast');
                %title('Joint Quantities','FontSize',15);
                set(leg,'Interpreter','latex');
                set(leg,'FontSize',13);
                xlabel('Time [s]','HorizontalAlignment','center',...
                       'FontWeight','bold',...
                       'FontSize',11,...
                       'Interpreter','latex');
                ylabel('$q$ [deg]','HorizontalAlignment','center',...
                       'FontWeight','bold',...
                       'FontSize',13,...
                       'Interpreter','latex');
                axis tight;
                grid on;  
                title('Joint quantities',...
                      'HorizontalAlignment','center',...
                      'Interpreter','latex');

                subplot(312);
                plot1 = plot(dataTime,(180/pi)*dq1,'lineWidth',1.5); hold on;
                set(plot1,'color',[1 0 0]);
                plot2= plot(dataTime,(180/pi)*dq2,'lineWidth',1.5); hold on;
                set(plot2,'color',[0 0.498039215803146 0]);
                leg = legend('$\dot q_{1}$','$\dot q_{2}$','Location','northeast');
                set(leg,'Interpreter','latex');
                set(leg,'FontSize',13);
                xlabel('Time [s]','HorizontalAlignment','center',...
                       'FontWeight','bold',...
                       'FontSize',11,...
                       'Interpreter','latex');
                ylabel('$\dot q$ [deg/s]','HorizontalAlignment','center',...
                       'FontWeight','bold',...
                       'FontSize',13,...
                       'Interpreter','latex');
                axis tight;
                grid on;

                subplot(313);
                plot1 = plot(dataTime,tau_ankle,'lineWidth',1.5); hold on;
                set(plot1,'color',[1 0 0]);
                plot2= plot(dataTime,tau_hip,'lineWidth',1.5); hold on;
                set(plot2,'color',[0 0.498039215803146 0]);
                leg = legend('$\tau_{1}$','$\tau_{2}$','Location','northeast');
                set(leg,'Interpreter','latex');
                set(leg,'FontSize',13);
                xlabel('Time [s]','HorizontalAlignment','center',...
                       'FontWeight','bold',...
                       'FontSize',11,...
                       'Interpreter','latex');
                ylabel('$\tau$ [Nm]','HorizontalAlignment','center',...
                       'FontWeight','bold',...
                       'FontSize',13,...
                       'Interpreter','latex');
                axis tight;
                grid on;
                
                if(plotPaper2)
                   save2pdf(fullfile(figFolder, ('JointQuantities')),fig,600);
                end
            end   
            fprintf('\n');
        end   
end 
   
%% ===============================SECTION3===============================%%
if(Section3)
        %% load
        load('./experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat');   
        load('./experiments/humanFixedBase/intermediateDataFiles/finalResults.mat');

        subjectList = 1;
        trialList = 2;
        
        for subjectID = subjectList
            fprintf('\n---------\nSubject : %d ',subjectID);
            for trialID = trialList
                fprintf('\nTrial : %d ',trialID);
                        
                currentTrial = finalResults(subjectID,trialID); 
                resMAP_2sens = currentTrial.resMAP_2sens;
                resMAP_3sens = currentTrial.resMAP_3sens;
                resMAP_4sens = currentTrial.resMAP_4sens;
                dataTime =currentTrial.dataTime;
                sample = 200;      %chose a generic sample to plot 
                
                if(plotMAP2torques)    
                %% Comparing MAP torques  
                    fig = figure();
                    axes1 = axes('Parent',fig,'FontSize',16);
                    box(axes1,'on');
                    hold(axes1,'on');
                    grid on;
                    %title(sprintf('Subject: %d, Trial: %d',subjectID, trialID))

                    plot1 = plot(dataTime,resMAP_4sens.tau_ankle, 'lineWidth',2.0,'LineStyle','-'); hold on;
                    set(plot1,'color',[1 0 0]);
                    plot2 = plot(dataTime,resMAP_4sens.tau_hip, 'lineWidth',2.0,'LineStyle','-'); hold on;
                    set(plot2,'color',[0 0.498039215803146 0]);

                    leg = legend('$\tau_{1}$','$\tau_{2}$','Location','southeast');
                    set(leg,'Interpreter','latex');
                    set(leg,'FontSize',18);
                    xlabel('Time [s]','FontSize',18);
                    ylabel('Joint torque [Nm]','FontSize',18);
                    axis tight;
                    grid on;
                

                %% Comparing MAP torques 1 and 2 in separate plots with associated variances (2sigma)

                    %====== TAU1
                    fig1 = figure();
                    axes1 = axes('Parent',fig1,'FontSize',16);
                    box(axes1,'on');
                    hold(axes1,'on');
                    grid on;
                    %title(sprintf('Subject: %d, Trial: %d',subjectID, trialID))
                    lineProps = {'LineWidth', 4.0};

                    %MAP 2 sens
                    lineProps.col = [0.87058824300766 0.490196079015732 0];
                    shad1 = shadedErrorBar(dataTime,resMAP_2sens.tau_ankle,2.*sqrt(resMAP_2sens.Stau_ankle),lineProps,1); 

                    %MAP 3 sens
                    lineProps.col = [0.494117647409439 0.184313729405403 0.556862771511078];
                    shad2 = shadedErrorBar(dataTime,resMAP_3sens.tau_ankle,2.*sqrt(resMAP_3sens.Stau_ankle),lineProps,1); 

                    %MAP 4 sens
                    lineProps.col = [0.466666668653488 0.674509823322296 0.18823529779911];
                    shad3 = shadedErrorBar(dataTime,resMAP_4sens.tau_ankle,2.*sqrt(resMAP_4sens.Stau_ankle),lineProps,1); 

                    leg = legend([shad1.mainLine,shad1.patch,shad2.mainLine,shad2.patch,shad3.mainLine,shad3.patch],...
                        {'$\mu_{\tau_1|{\ddot q}}$','$2\sigma_{\tau_1|{\ddot q}}$','$\mu_{\tau_1|{\ddot q,y_1}}$','$2\sigma_{\tau_1|{\ddot q,y_1}}$','$\mu_{\tau_1|{\ddot q,y_1,y_2}}$','$2\sigma_{\tau_1|{\ddot q,y_1,y_2}}$'},'Location','northeast');
                    set(leg,'Interpreter','latex');
                    set(leg,'FontSize',18);

                    xlabel('Time [s]','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',18,...
                           'Interpreter','latex');
                    ylabel('$\tau_1$ [Nm]','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',18,...
                           'Interpreter','latex');
                    axis tight;
                    box(axes1,'on');
                    grid on;

                    if(plotPaperTau1)
                        save2pdf(fullfile(figFolder, ('tau1')),fig1,600);
                    end

                    %====== TAU2
                    fig2 = figure();
                    axes1 = axes('Parent',fig2,'FontSize',16);
                    box(axes1,'on');
                    hold(axes1,'on');
                    grid on;
                    %title(sprintf('Subject: %d, Trial: %d',subjectID, trialID))
                    lineProps = {'LineWidth', 4.0};

                    %MAP 2 sens
                    lineProps.col = [0.87058824300766 0.490196079015732 0];
                    shad1 = shadedErrorBar(dataTime,resMAP_2sens.tau_hip,2.*sqrt(resMAP_2sens.Stau_hip),lineProps,1);

                    %MAP 3 sens
                    lineProps.col = [0.494117647409439 0.184313729405403 0.556862771511078];
                    shad2 = shadedErrorBar(dataTime,resMAP_3sens.tau_hip,2.*sqrt(resMAP_3sens.Stau_hip),lineProps,1);

                    %MAP 4 sens
                    lineProps.col = [0.466666668653488 0.674509823322296 0.18823529779911];
                    shad3 = shadedErrorBar(dataTime,resMAP_4sens.tau_hip,2.*sqrt(resMAP_4sens.Stau_hip),lineProps,1);

                    leg = legend([shad1.mainLine,shad1.patch,shad2.mainLine,shad2.patch,shad3.mainLine,shad3.patch],...
                        {'$\mu_{\tau_2|{\ddot q}}$','$2\sigma_{\tau_2|{\ddot q}}$','$\mu_{\tau_2|{\ddot q,y_1}}$','$2\sigma_{\tau_2|{\ddot q,y_1}}$','$\mu_{\tau_2|{\ddot q,y_1,y_2}}$','$2\sigma_{\tau_2|{\ddot q,y_1,y_2}}$'},'Location','northeast');
                    set(leg,'Interpreter','latex');
                    set(leg,'FontSize',18);

                    xlabel('Time [s]','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',18,...
                           'Interpreter','latex');
                    ylabel('$\tau_2$ [Nm]','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',18,...
                           'Interpreter','latex');
                    axis tight;
                    grid on;

                    if(plotPaperTau2)
                        save2pdf(fullfile(figFolder, ('tau2')),fig2,600);
                    end 
                end 
                
                %% Covariances (complete Sd 52x52)
                if(plotCovariances)
     
                    clims = [-4,6.5];  %limits of the colorbar
                    labels = {'a_1','fB_1','f_1','\tau_1','fx_1','ddq_1','a_2','fB_2','f_2','\tau_2','fx_2','ddq_2'};

                    fig1 = figure();
                    axes1 = axes('Parent',fig1,...
                                 'XTickLabel',labels,...          % x labels
                                 'XTick',resMAP_2sens.jIndex ,...  % x indices
                                 'XTickLabelRotation',90,...      % x 90deg labels rotation
                                 'XAxisLocation','top',...        % x location 
                                 'YTickLabel',labels,...          % y labels
                                 'YTick',resMAP_2sens.jIndex ,...  % y indices
                                 'Layer','top',...
                                 'YDir','reverse',...
                                 'FontSize',12,...                % fontsize of ticks
                                 'DataAspectRatio',[1 1 1]);      % proportional tration of the image
                    hold(axes1,'on');
                    imagesc(resMAP_2sens.Sd{sample,1},clims);
                    %title('MAP 2 sensors (ddq & ftx)', 'FontSize',15);
                    axis tight;
                    box(axes1,'on');
                    colorbar;
                    if(plotPaperCov)
                        save2pdf(fullfile(figFolder, ('covariance2sens')),fig1,600);
                    end


                    fig2 = figure();
                    axes1 = axes('Parent',fig2,...
                                 'XTickLabel',labels,...          % x labels
                                 'XTick',resMAP_2sens.jIndex ,...  % x indices
                                 'XTickLabelRotation',90,...      % x 90deg labels rotation
                                 'XAxisLocation','top',...        % x location 
                                 'YTickLabel',labels,...          % y labels
                                 'YTick',resMAP_2sens.jIndex ,...  % y indices
                                 'Layer','top',...
                                 'YDir','reverse',...
                                 'FontSize',12,...                % fontsize of ticks
                                 'DataAspectRatio',[1 1 1]);      % proportional tration of the image
                    hold(axes1,'on');
                    imagesc(resMAP_3sens.Sd{sample,1},clims);
                    %title('MAP 3 sensors (ddq & ftx & fp)', 'FontSize',15);
                    axis tight;
                    box(axes1,'on');
                    colorbar;
                    if(plotPaperCov)
                        save2pdf(fullfile(figFolder, ('covariance3sens')),fig2,600);
                    end


                    fig3 = figure();
                    axes1 = axes('Parent',fig3,...
                                 'XTickLabel',labels,...          % x labels
                                 'XTick',resMAP_2sens.jIndex ,...  % x indices
                                 'XTickLabelRotation',90,...      % x 90deg labels rotation
                                 'XAxisLocation','top',...        % x location 
                                 'YTickLabel',labels,...          % y labels
                                 'YTick',resMAP_2sens.jIndex ,...  % y indices
                                 'Layer','top',...
                                 'YDir','reverse',...
                                 'FontSize',12,...                % fontsize of ticks
                                 'DataAspectRatio',[1 1 1]);      % proportional tration of the image
                    hold(axes1,'on');
                    imagesc(resMAP_4sens.Sd{sample,1},clims);
                    %title('MAP all sensors)', 'FontSize',15);
                    axis tight;
                    box(axes1,'on');
                    colorbar;
                    if(plotPaperCov)
                        save2pdf(fullfile(figFolder, ('covariance4sens')),fig3,600);
                    end
                end
            
                %% Bar Covariances
                if(plotBarCovariances) 
              
%                 % =====d2q
%                 fig1 = figure();
%                 axes1 = axes('Parent',fig1,'FontSize',16);
%                 box(axes1,'on');
%                 hold(axes1,'on');
%                 grid on;
%                 title('d2q')
%                 
%                 Sd_d2q_2sens = diag ([resMAP_2sens.Sd2q_ankle(sample),resMAP_2sens.Sd2q_hip(sample)]);
%                 Sd_d2q_3sens = diag ([resMAP_3sens.Sd2q_ankle(sample),resMAP_3sens.Sd2q_hip(sample)]);
%                 Sd_d2q_4sens = diag ([resMAP_4sens.Sd2q_ankle(sample),resMAP_4sens.Sd2q_hip(sample)]);
%   
%                 Ad2q1 = diag(Sd_d2q_2sens);
%                 Ad2q2 = diag(Sd_d2q_3sens);
%                 Ad2q3 = diag(Sd_d2q_4sens);
%                 
% %                 plot(Ad2q1,'b');hold on;
% %                 plot(Ad2q2,'r');
% %                 plot(Ad2q3,'g');
%                 
% %                 figure();
%                 bar([1, 2],[Ad2q1(1) Ad2q1(2)],0.3,'b'); hold on;
%                 bar([1.3, 2.3],[Ad2q2(1) Ad2q2(2)],0.3,'r');
%                 bar([1.6, 2.6],[Ad2q3(1) Ad2q3(2)],0.3,'g');
% 
%                 leg = legend('$\d2q_{2sens}$','$\d2q_{3sens}$','$\d2q_{4sens}$','Location','southeast');
%                 set(leg,'Interpreter','latex');
%                 set(leg,'FontSize',18);
%                 xlabel('xxxx','FontSize',20);
%                 ylabel('xxxx','FontSize',20);
%                 %axis tight;
%                 grid on;
                
                
                    % =====TAU
                    fig1 = figure();
                    axes1 = axes('Parent',fig1,...
                                 'YGrid','on',...
                                 'YTick',[0 1 2 3 4],...
                                 'XTickLabel',{'$\tau_{1}$','$\tau_{2}$'},...
                                 'XTick',[1.3 2.3 ],...
                                 'TickLabelInterpreter','latex',...
                                 'FontSize',18);
                    box(axes1,'off');
                    hold(axes1,'on');

    %                 title('Joint torques covariances',...
    %                       'HorizontalAlignment','center',...
    %                       'Interpreter','latex');

                    Sd_tau_2sens = diag ([sqrt(resMAP_2sens.Stau_ankle(sample)),sqrt(resMAP_2sens.Stau_hip(sample))]);
                    Sd_tau_3sens = diag ([sqrt(resMAP_3sens.Stau_ankle(sample)),sqrt(resMAP_3sens.Stau_hip(sample))]);
                    Sd_tau_4sens = diag ([sqrt(resMAP_4sens.Stau_ankle(sample)),sqrt(resMAP_4sens.Stau_hip(sample))]);

                    Atau1 = diag(Sd_tau_2sens);
                    Atau2 = diag(Sd_tau_3sens);
                    Atau3 = diag(Sd_tau_4sens);

                    bar1 = bar([1, 2],[Atau1(1) Atau1(2)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 = bar([1.3, 2.3],[Atau2(1) Atau2(2)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 = bar([1.6, 2.6],[Atau3(1) Atau3(2)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]);

                    leg = legend('$\ddot q$','$\ddot q,y_1$','$\ddot q,y_1,y_2$','Location','east');
                    set(leg,'Interpreter','latex');
                    set(leg,'FontSize',18);                
                    xlabel('joint torques','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',20,...
                           'Interpreter','latex');
                    ylabel('standard deviation [Nm]','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',20,...
                           'Interpreter','latex');

                     if(plotPaperBarCov)
                        save2pdf(fullfile(figFolder, ('barCovariancesTau')),fig1,600);
                     end


                    % =====MOMENTS
                    fig2 = figure();
                    title('Joint torques covariances');

                    Sf_leg = diag(sqrt(resMAP_2sens.Sf_leg{sample}));
                    Sf_torso = diag(sqrt(resMAP_2sens.Sf_torso{sample}));
                    Sd_force_2sens = diag ([Sf_leg',Sf_torso']);
                    Aforce1 = diag(Sd_force_2sens);

                    Sf_leg = diag(sqrt(resMAP_3sens.Sf_leg{sample}));
                    Sf_torso = diag(sqrt(resMAP_3sens.Sf_torso{sample}));
                    Sd_force_3sens = diag ([ Sf_leg',Sf_torso']);
                    Aforce2 = diag(Sd_force_3sens);

                    Sf_leg = diag(sqrt(resMAP_4sens.Sf_leg{sample}));
                    Sf_torso = diag(sqrt(resMAP_4sens.Sf_torso{sample}));
                    Sd_force_4sens = diag ([Sf_leg',Sf_torso']);
                    Aforce3 = diag(Sd_force_4sens);


                    subplot1 = subplot(3,1,1,'Parent',fig2,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 2.5 5],...
                               'XTickLabel',{'$\mu_{1}(1)$','$\mu_{2}(1)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot1,'off');
                    hold(subplot1,'on');
                    bar1 = bar([1, 2],[Aforce1(1) Aforce1(7)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 =bar([1.3, 2.3],[Aforce2(1) Aforce2(7)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 =bar([1.6, 2.6],[Aforce3(1) Aforce3(7)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]);                
                    ylim(subplot1,[0 3]);

                    
                    subplot2 = subplot(3,1,2,'Parent',fig2,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 2.5 5],...
                               'XTickLabel',{'$\mu_{1}(2)$','$\mu_{2}(2)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot2,'off');
                    hold(subplot2,'on');
                    bar1 = bar([1, 2],[Aforce1(2) Aforce1(8)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 = bar([1.3, 2.3],[Aforce2(2) Aforce2(8)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 = bar([1.6, 2.6],[Aforce3(2) Aforce3(8)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]); 
                    ylim(subplot2,[0 3]);
                    ylabel('standard deviation [Nm]','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',20,...
                           'Interpreter','latex');


                    subplot3 = subplot(3,1,3,'Parent',fig2,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 2.5 5],...
                               'XTickLabel',{'$\mu_{1}(3)$','$\mu_{2}(3)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot3,'off');
                    hold(subplot3,'on');
                    bar1 = bar([1, 2],[Aforce1(3) Aforce1(9)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 = bar([1.3, 2.3],[Aforce2(3) Aforce2(9)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 = bar([1.6, 2.6],[Aforce3(3) Aforce3(9)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]); 
                    ylim(subplot3,[0 3]);
                    xlabel('couple components','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',20,...
                           'Interpreter','latex');

                    if(plotPaperBarCov)
                        save2pdf(fullfile(figFolder, ('barCovariancesMoment')),fig2,600);
                    end

                    
                    % =====LINEAR FORCE
                    fig3 =figure();

                    subplot1 = subplot(3,1,1,'Parent',fig3,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 2.5 5],...
                               'XTickLabel',{'$f_{1}(1)$','$f_{2}(1)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot1,'off');
                    hold(subplot1,'on');
                    bar1 = bar([1, 2],[Aforce1(4) Aforce1(10)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 =bar([1.3, 2.3],[Aforce2(4) Aforce2(10)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 =bar([1.6, 2.6],[Aforce3(4) Aforce3(10)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]);
                    ylim(subplot1,[0 3]);

                    
                    subplot2 = subplot(3,1,2,'Parent',fig3,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 2.5 5],...
                               'XTickLabel',{'$f_{1}(2)$','$f_{2}(2)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot2,'off');
                    hold(subplot2,'on');
                    bar1 = bar([1, 2],[Aforce1(5) Aforce1(11)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 = bar([1.3, 2.3],[Aforce2(5) Aforce2(11)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 = bar([1.6, 2.6],[Aforce3(5) Aforce3(11)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]); 
                    ylim(subplot2,[0 3]);
                    ylabel('standard deviation [N] ','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',20,...
                           'Interpreter','latex');


                    subplot3 = subplot(3,1,3,'Parent',fig3,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 2.5 5],...
                               'XTickLabel',{'$f_{1}(3)$','$f_{2}(3)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot3,'off');
                    hold(subplot3,'on');
                    bar1 = bar([1, 2],[Aforce1(6) Aforce1(12)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 = bar([1.3, 2.3],[Aforce2(6) Aforce2(12)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 = bar([1.6, 2.6],[Aforce3(6) Aforce3(12)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]); 
                    ylim(subplot3,[0 3]);
                    xlabel('linear force components','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',20,...
                           'Interpreter','latex');

                    if(plotPaperBarCov)
                        save2pdf(fullfile(figFolder, ('barCovariancesForce')),fig3,600);
                    end    

                    
                    % =====ANGULAR ACC
                    fig4 = figure();
                    title('Joint torques covariances');

                    Sa_leg = diag(sqrt(resMAP_2sens.Sa_leg{sample}));
                    Sa_torso = diag(sqrt(resMAP_2sens.Sa_torso{sample}));
                    Sd_acc_2sens = diag ([ Sa_leg',Sa_torso']);
                    Aacc1 = diag(Sd_acc_2sens);

                    Sa_leg = diag(sqrt(resMAP_3sens.Sa_leg{sample}));
                    Sa_torso = diag(sqrt(resMAP_3sens.Sa_torso{sample}));
                    Sd_acc_3sens = diag ([  Sa_leg',Sa_torso']);
                    Aacc2 = diag(Sd_acc_3sens);

                    Sa_leg = diag(sqrt(resMAP_4sens.Sa_leg{sample}));
                    Sa_torso = diag(sqrt(resMAP_4sens.Sa_torso{sample}));
                    Sd_acc_4sens = diag ([  Sa_leg',Sa_torso']);
                    Aacc3 = diag(Sd_acc_4sens);


                    subplot1 = subplot(3,1,1,'Parent',fig4,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 0.05],...
                               'XTickLabel',{'$\dot \omega_{1}(1)$','$\dot \omega_{2}(1)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot1,'off');
                    hold(subplot1,'on');
                    bar1 = bar([1, 2],[Aacc1(1) Aacc1(7)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 =bar([1.3, 2.3],[Aacc2(1) Aacc2(7)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 =bar([1.6, 2.6],[Aacc3(1) Aacc3(7)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]);                
                    ylim(subplot1,[0 0.05]);

                    
                    subplot2 = subplot(3,1,2,'Parent',fig4,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 0.05],...
                               'XTickLabel',{'$\dot \omega_{1}(2)$','$\dot \omega_{2}(2)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot2,'off');
                    hold(subplot2,'on');
                    bar1 = bar([1, 2],[Aacc1(2) Aacc1(8)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 = bar([1.3, 2.3],[Aacc2(2) Aacc2(8)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 = bar([1.6, 2.6],[Aacc3(2) Aacc3(8)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]); 
                    ylim(subplot2,[0 0.05]);
                    ylabel('standard deviation [$rad/s^2$]','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',20,...
                           'Interpreter','latex');


                    subplot3 = subplot(3,1,3,'Parent',fig4,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 0.05],...
                               'XTickLabel',{'$\dot \omega_{1}(3)$','$\dot \omega_{2}(3)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot3,'off');
                    hold(subplot3,'on');
                    bar1 = bar([1, 2],[Aacc1(3) Aacc1(9)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 = bar([1.3, 2.3],[Aacc2(3) Aacc2(9)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 = bar([1.6, 2.6],[Aacc3(3) Aacc3(9)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]); 
                    ylim(subplot3,[0 0.05]);
                    xlabel('angular acceleration components','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',20,...
                           'Interpreter','latex');

                    if(plotPaperBarCov)
                        save2pdf(fullfile(figFolder, ('barCovariancesAngAcc')),fig4,600);
                    end

                    
                    % =====LINEAR ACC
                    fig5 =figure();

                    subplot1 = subplot(3,1,1,'Parent',fig5,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 0.05],...
                               'XTickLabel',{'$a_{1}(1)$','$a_{2}(1)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot1,'off');
                    hold(subplot1,'on');
                    bar1 = bar([1, 2],[Aacc1(4) Aacc1(10)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 =bar([1.3, 2.3],[Aacc2(4) Aacc2(10)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 =bar([1.6, 2.6],[Aacc3(4) Aacc3(10)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]);
                    ylim(subplot1,[0 0.05]);

                    
                    subplot2 = subplot(3,1,2,'Parent',fig5,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 0.05],...
                               'XTickLabel',{'$a_{1}(2)$','$a_{2}(2)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot2,'off');
                    hold(subplot2,'on');
                    bar1 = bar([1, 2],[Aacc1(5) Aacc1(11)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 = bar([1.3, 2.3],[Aacc2(5) Aacc2(11)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 = bar([1.6, 2.6],[Aacc3(5) Aacc3(11)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]); 
                    ylim(subplot2,[0 0.05]);
                    ylabel('standard deviation [$m/s^2$] ','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',20,...
                           'Interpreter','latex');


                    subplot3 = subplot(3,1,3,'Parent',fig5,'YGrid','on','XGrid','off',...
                               'TitleFontSizeMultiplier',10,...
                               'YTick',[0 0.05],...
                               'XTickLabel',{'$a_{1}(3)$','$a_{2}(3)$'},...
                               'XTick',[1.3 2.3],...
                               'TickLabelInterpreter','latex',...
                               'FontWeight','bold',...
                               'FontSize',18);
                    box(subplot3,'off');
                    hold(subplot3,'on');
                    bar1 = bar([1, 2],[Aacc1(6) Aacc1(12)],0.3,'FaceColor',[0.87058824300766 0.490196079015732 0]); hold on;
                    bar2 = bar([1.3, 2.3],[Aacc2(6) Aacc2(12)],0.3,'FaceColor',[0.494117647409439 0.184313729405403 0.556862771511078]);
                    bar3 = bar([1.6, 2.6],[Aacc3(6) Aacc3(12)],0.3,'FaceColor',[0.466666668653488 0.674509823322296 0.18823529779911]); 
                    ylim(subplot3,[0 0.05]);
                    xlabel('linear acceleration components','HorizontalAlignment','center',...
                           'FontWeight','bold',...
                           'FontSize',20,...
                           'Interpreter','latex');

                    if(plotPaperBarCov)
                        save2pdf(fullfile(figFolder, ('barCovariancesLinAcc')),fig5,600);
                    end    
                end   
            end   
        end 
end

%% ===============================SECTION4===============================%%
if(Section4)
        %% load data
        load('./experiments/humanFixedBase/intermediateDataFiles/synchronisedData.mat');   

        subjectList = 1;
        trialList = 3;
        
        for subjectID = subjectList
            fprintf('\n---------\nSubject : %d ',subjectID);
            for trialID = trialList
                fprintf('\nTrial : %d ',trialID);

                currentSynchTrial = synchronisedData(subjectID,trialID);
                dataTime = currentSynchTrial.dataTime;
                aLin_imu_imu = currentSynchTrial.aLin_imu_imu;
                omega_imu_imu = currentSynchTrial.omega_imu_imu;
                wrench_fp_fp = currentSynchTrial.wrench_fp_fp;

                %% IMU data (a_imuLin, omega_imu)

                % Plotting raw data coming from IMU sensor in IMU frame     
                fig1 = figure();
                axes1 = axes('Parent',fig1,'FontSize',16);
                box(axes1,'on');
                hold(axes1,'on');
                grid on; 

                subplot(211);
                plot(dataTime, aLin_imu_imu); axis tight;
                leg = legend('$a_x$','$a_y$','$a_z$','Location','northeast');
                set(leg,'Interpreter','latex');
                set(leg,'FontSize',15);
                xlabel('Time [s]','FontSize',15);
                ylabel('Linear Acceleration [m/sec^2]','FontSize',15);
                %title('Raw IMU data of link 2 (IMU frame)','FontSize',15);
                grid on;
                title(sprintf('Subject %d, Trial %d, Raw IMU data of link 2 (IMU frame)',subjectID,trialID));

                subplot(212);
                plot(dataTime,omega_imu_imu);  axis tight;
                leg = legend('$w_x$','$w_y$','$w_z$','Location','northeast');
                set(leg,'Interpreter','latex');
                set(leg,'FontSize',15);
                xlabel('Time [s]','FontSize',15);
                ylabel('Angular Velocity [rad/s]','FontSize',15);
                grid on;
                
                if(plotSensor)
                        save2pdf(fullfile(figFolder, ('IMUsensor')),fig1,600);
                end

                %% forceplate wrench data (f_fp_fp)

                % Plotting raw data coming from force plate sensor in sensor frame      
                fig2 = figure();
                axes1 = axes('Parent',fig2,'FontSize',16);
                box(axes1,'on');
                hold(axes1,'on');
                grid on;

                subplot(211);
                plot(dataTime,wrench_fp_fp(:,4:6)); axis tight;
                leg =legend('$F_x$','$F_y$','$F_z$','Location','northeast');
                set(leg,'Interpreter','latex');
                set(leg,'FontSize',15);
                xlabel('Time [s]','FontSize',15);
                ylabel('Force [N]','Fontsize',15);
                %title('Wrench measured in force plate (Fp frame)','FontSize',15);
                grid on;
                title(sprintf('Subject %d, Trial %d, Wrench measured in force plate (Fp frame)',subjectID,trialID));

                subplot(212);
                plot(dataTime,wrench_fp_fp(:,1:3)); axis tight;
                leg = legend('$M_x$','$M_y$','$M_z$','Location','northeast');
                set(leg,'Interpreter','latex');
                set(leg,'FontSize',15);
                xlabel('Time [s]','FontSize',15);
                ylabel('Moment [Nm]','FontSize',15);
                grid on; 
                
                if(plotSensor)
                        save2pdf(fullfile(figFolder, ('ForcePlateSensor')),fig2,600);
                end
            end
        end
end