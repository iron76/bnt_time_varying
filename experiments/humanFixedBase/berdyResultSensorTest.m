%% load dataset
subjectID = 1; trialID = 1;
load(sprintf('./experiments/humanFixedBase/savedBERDYresult_subj%d_trial%d.mat',subjectID,trialID),'res','data','myMAP');

y_predRNEA = zeros(26,length(data.time));
q = zeros(2,(length(data.time)));dq = zeros(2,(length(data.time)));ddq = zeros(2,(length(data.time)));
q(1,:) = data.q1;q(2,:) = data.q2;
dq(1,:) = data.dq1;dq(2,:) = data.dq2;
ddq(1,:) = data.ddq1;ddq(2,:) = data.ddq2;

load(sprintf('./experiments/humanFixedBase/humanThreeLinkModelFromURDF_subject%d.mat',subjectID));
dmodel  = humanThreeLink_dmodel; 

fprintf('\nComputing RNEA prediction\n');
%% time loop to extract RNEA simY
for tIdx = 1:length(data.time)
    f_ext={zeros(6,1),zeros(6,1)};
    y_predRNEA(:,tIdx) = simulateSensors_RNEA(dmodel,q(:,tIdx),dq(:,tIdx),ddq(:,tIdx),f_ext);
end


%% indices in output chosen for comparison and analysis
chosenInd = [1:12,25,26];
% [f mu a2lin a2rot fx1 fx2 dq1 dq2]
% [3  3  3      3    6   6   1   1 ]

fprintf('\nComputing MAP prediction\n');
%% test Y from MAP simY
y_predMAP = myMAP.simY(res.d);

fprintf('\nPlotting selected data\n');
%% chosen output indices
 for  ind = chosenInd
     
        if(ind>=1 && ind<7)
            titl = 'fts';
        elseif(ind>6 && ind<13)
            titl = 'imu';
        elseif(ind>12 && ind<25)
            titl = 'fx';
        elseif(ind>24 && ind<27)
            titl = 'ddq';
        end

        fig = figure();
        axes1 = axes('Parent',fig,'FontSize',16);
        box(axes1,'on');
        hold(axes1,'on');
        grid on;

        plot1 = plot(data.time,y_predMAP(ind,:), 'lineWidth',1.0, 'LineStyle','--'); hold on;
        set(plot1,'color',[1 0 0]);
        plot2 = plot(data.time,y_predRNEA(ind,:), 'lineWidth',1.0, 'LineStyle','-.'); hold on;
        set(plot2,'color',[0 1 0]);
        plot3 = plot(data.time,data.y(ind,:), 'lineWidth',1.0); hold on;
        set(plot3,'color',[0 0 1]);

        leg = legend('MAP Pred', 'RNEA Pred' , 'Actual data','Location','northeast');
        %set(leg,'Interpreter','latex');
        set(leg,'FontSize',18);
        xlabel('Time [s]','FontSize',20);
        %ylabel('Torque[Nm]','FontSize',20);
        title(sprintf('Y[%d] : %s',ind,titl));
        axis tight;
        grid on;
 end
