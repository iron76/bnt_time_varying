clear;clc;close all;

%% testOptions
plotMAPtorques = true;
plotCovariances = false;
plotPaper = false;

%% folder for plots
figFolder = './experiments/humanFixedBase/plots';
if(exist(figFolder,'dir')==0)
    mkdir(figFolder);
end

%% selected subjects and trials
subjectList = 8;
trialList = 1;  

for subjectID = subjectList
    fprintf('\n---------\nSubject : %d ',subjectID);
    for trialID = trialList
        fprintf('\nTrial : %d\n',trialID);
        
        %% load BERDY compatible data and sensor link transforms
        load('./experiments/humanFixedBase/intermediateDataFiles/BERDYFormattedSensorData.mat');

        currentTrial = BERDYFormattedSensorData(subjectID,trialID); 
        data = currentTrial.data;
        dataTime = currentTrial.data.dataTime;
        q = currentTrial.data.q';
        dq = currentTrial.data.dq';
        ddq = currentTrial.data.ddq';
        len = length(dataTime);

        load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');

        currentTrialSens = sensorLinkTransforms(subjectID,trialID);
        X_imu_2 = currentTrialSens.X_imu_2;
        XStar_fp_0 = currentTrialSens.XStar_fp_0;
        XStar_0_1 = currentTrialSens.XStar_0_1;

        % =====structure from files
        data.parts    = {'leg'         ,'torso'};
        data.labels   = {'fts'         ,'imu'  };
        data.ndof     = {6             ,6      };
        data.index    = {'1:6'         ,'1:6'  };

        % =====structure of sensors
        sens.parts    = {'leg'         ,'torso'};       %force of the forceplate is trasmitted into the leg
        sens.labels   = {'fts'         ,'imu'  };
        sens.ndof     = {6             ,6      };

        label_to_plot = {'fts'         ,'imu'  };

        %% build models
        load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');

        currentModel = humanThreeLinkModelFromURDF(subjectID);

        humanThreeLink_dmodel.linkname = {'leg' 'torso'}; 
        humanThreeLink_dmodel.jointname = {'ankle' 'hip'}; 

        dmodel  = currentModel.dmodel;                      %deterministic model 
        ymodel  = humanThreeLinkSens(dmodel, sens);  

        dmodel  = autoTreeStochastic(dmodel, 1e-6);         % probabilistic model for D equation (added Sv and Sw)
        ymodel  = humanThreeLinkSensStochastic(ymodel);     % probabilistic model for Y(q,dq) d = y (added Sy)

        myModel = model(dmodel);
        mySens  = sensors(ymodel);  

        myMAP   = MAP(myModel, mySens);
   
        %% ================================ RNEA ==================================

        ymodel_RNEA  = autoSensRNEA(dmodel);
        mySens_RNEA  = sensors(ymodel_RNEA);
        myRNEA       = RNEA(myModel, mySens_RNEA);

        y_RNEA_f = zeros(6*dmodel.NB, len);
        y_RNEA_ddq = zeros(dmodel.NB, len);
        fx = cell(dmodel.NB);

        %Ordering y_RNEA in the form [fx1 fx2 ddq1 ddq2]
        for i = 1 : dmodel.NB
            for t = 1 : len
                fx{i,1} = zeros(6,1); 
                y_RNEA_f(((1:6)+(i-1)*6), t) = [fx{i,1}];
                y_RNEA_ddq(i, t) = [ddq(i,t)];
            end
            y_RNEA = [y_RNEA_f ; y_RNEA_ddq];
        end


        d= zeros (26*myRNEA.IDmodel.modelParams.NB,len);
        for i = 1 : len
             myRNEA = myRNEA.setState(q(:,i), dq(:,i));
             myRNEA = myRNEA.setY(y_RNEA(:,i));
             myRNEA = myRNEA.solveID();

             d(:,i) = myRNEA.d;
             
             resRNEA.tau(i,:) = myRNEA.tau;
        end
        
        clear fx;
        % ========end RNEA
        
        %% Build data.y anda data.Sy 
        % data.y are ordered in:  
        % - angular-linear notation 
        % - the form [f1 a2 ftx1 ftx2 ddq1 ddq2]
        % - sensor frame

        %===== data.y
        data.y  = [];
        for i = 1 : length(sens.labels)
             eval(['data.y  = [data.y  data.ys_sensFrame_' sens.labels{i} '];']);
        end
        % Add null external forces ftx = 0
        data.y  = [data.y, zeros(len,6*dmodel.NB)];
        % Add ddq measurements
        data.y  = [data.y, data.ddq];
        data.y = data.y';

        %===== data.Sy
        data.Sy = [];
        for i = 1 : length(myMAP.IDsens.sensorsParams.labels)
             data.Sy = [data.Sy; diag(myMAP.IDsens.sensorsParams.Sy{i})];
        end
        data.Sy = repmat(data.Sy, 1, len-1);
        data.Sy = [data.Sy data.Sy(:,end)];

        %% ================================ MAP ===================================
        %% Build Ymatrix manually
        % Ymatrix has to be consistent with measurements form [f1 a2 ftx1 ftx2 ddq1 ddq2]

        Y = cell2mat(ymodel.Y);
        %Y = zeros (ymodel.m,26*dmodel.NB);
        Y(7:12,27:32) = X_imu_2;
        Y(13:18,20:25) = eye(6);
        Y(19:24,46:51) = eye(6);
        Y(25,26) = eye(1);
        Y(26,52) = eye(1);


        Ymatrix = cell(len,1);
        for i = 1 : len
            %the only row in Ymatrix that is time varying
            Y(1:6,13:18) = XStar_fp_0 * XStar_0_1{i};
            Ymatrix{i} = Y; 
        end

        clear Y;
   
        %% Build bias b_Y manually
        % b_Y has to be consistent with Ymatrix

        load('./experiments/humanFixedBase/data/subjectSizeParams.mat');

        currentParams = subjectParams(subjectID);

        footMass =  currentParams.footMass;
        posP_0 = [0; 0; (0.5*currentParams.footHeight)];
        footIxx =  currentParams.footIxx;
        footIyy =  currentParams.footIyy;
        footIzz =  currentParams.footIzz;

        b_Y = zeros (size(data.y)); 
        R_imu_2 = X_imu_2(1:3,1:3);

        a_G_grav = [0;0;0;0;0;-9.8100]; %Featherstone-like notation, in global reference
        X_G_0 = currentTrialSens.X_G_0;
        X_0_G = InverseAdjTransform(X_G_0);
        a_0_grav = X_0_G * a_G_grav;

        I_0 = createSpatialInertia(footIxx,footIyy,footIzz,footMass,posP_0);

        b_Y(1:6,1:len)   = repmat((-XStar_fp_0 * I_0 * a_0_grav),1,len);
        
        v = cell(size(q))';
        fx = zeros (6,1);
        fext    = cell(1,2);
        for i = 1 : dmodel.NB
             fext{i}    = fx;
        end

        
        % For the following computation are necessary velocities.  For this
        %reason it is utiliseda modified version of ID of Featherstone.  
        % So the path is not:
        % ../bnt_time_varying/extern/featherstone/dynamics/ID.m
        % but
        % ../bnt_time_varying/experiments/humanFixedBase/helperFunction/IDv.m
        for i = 1 : len 
            [~, ~,v_i,~,~] = IDv(dmodel, q(:,i), dq(:,i), ddq(:,i), fext); 
            v(i,:) = v_i;

            A =R_imu_2*v{i,2}(1:3,1);
            B =((X_imu_2(4:6,1:3)*v{i,2}(1:3,1))+(R_imu_2*v{i,2}(4:6,1)));
            b_Y(10:12,i) = cross(A,B);
        end 
        
        clear A;
        clear B;
        clear footIxx;
        clear footIyy;
        clear footIzz;
        
        %% Computing MAP method

        resMAP.Sd = cell(len,1);
        
        for i = 1 : len

            myMAP = myMAP.setState(data.q(i,:)', data.dq(i,:)');
            myMAP = myMAP.setY(data.y(:,i));
            myMAP = myMAP.setYmatrix(Ymatrix{i});
            myMAP = myMAP.setBias(b_Y(:,i));
            myMAP = myMAP.solveID();

           
            resMAP.d(:,i)       = myMAP.d;
            %resMAP.Sd(:,:,i)    = myMAP.Sd; %full() passing from sparse to double matrix
            resMAP.Sd{i} = myMAP.Sd;
            resMAP.Ymatrix{i,1} = myMAP.IDsens.sensorsParams.Y; 
            resMAP.b_Y(:,i)     = myMAP.IDsens.sensorsParams.bias;
            %resMAP.y(:,i)      = (Ymatrix{i} * res.d(:,i)) + b_Y(:,i); 
            %resMAP.Sy(:,:,i)   = Ymatrix{i} * res.Sd(:,:,i) * Ymatrix{i}';
            
             if mod(i-1,100) == 0
                    fprintf('Processing %d %% of the dataset\n', round(i/len*100));
             end
        end
        % ========end MAP
        %% Rerrange solution

        for i = 1 : dmodel.NB
             for j = 1 : len

                 link = strrep(myMAP.IDmodel.modelParams.linkname{i}, '+', '_');
                 joint = strrep(myMAP.IDmodel.modelParams.jointname{i}, '+', '_');

                 di   = ['resMAP.d_'   link  '(:,j)'];
                 ind  = '1 + 26*(i-1) : 26*(i-1) + 26';
                 eval([di '   = d(' ind ',j);'])

                 %a
                 ind  = '1 + 26*(i-1) : 26*(i-1) +  6';
                 eval(['resMAP.a_'    link '(:,j) = resMAP.d      (' ind '        ,j);'])
                 eval(['resMAP.Sa_'   link '{j,1} = resMAP.Sd{j,1}(' ind ',' ind '  );'])
                 %fB
                 ind  = '7 + 26*(i-1) : 26*(i-1) + 12';
                 eval(['resMAP.fB_'   link '(:,j) = resMAP.d      (' ind '        ,j);'])
                 eval(['resMAP.SfB_'  link '{j,1} = resMAP.Sd{j,1}(' ind ',' ind '  );'])
                 %f
                 ind  = '13 + 26*(i-1) : 26*(i-1) + 18';
                 eval(['resMAP.f_'    link '(:,j) = resMAP.d      (' ind '        ,j);'])
                 eval(['resMAP.Sf_'   link '{j,1} = resMAP.Sd{j,1}(' ind ',' ind '  );'])
                 %tau
                 ind  = '19 + 26*(i-1) : 26*(i-1) + 19';
                 eval(['resMAP.tau_'  joint '(:,j) = resMAP.d      (' ind '        ,j);'])
                 eval(['resMAP.Stau_' joint '{j,1} = resMAP.Sd{j,1}(' ind ',' ind '  );'])
                 %fx
                 ind  = '20 + 26*(i-1) : 26*(i-1) + 25';
                 eval(['resMAP.fx_'   link '(:,j) = resMAP.d      (' ind '        ,j);'])
                 eval(['resMAP.Sfx_'  link '{j,1} = resMAP.Sd{j,1}(' ind ',' ind '  );'])
                 %ddq
                 ind  = '26 + 26*(i-1) : 26*(i-1) + 26';
                 eval(['resMAP.d2q_'  joint '(:,j) = resMAP.d      (' ind '        ,j);'])
                 eval(['resMAP.Sd2q_' joint '{j,1} = resMAP.Sd{j,1}(' ind ',' ind '  );'])
             end
        end 
    
        %% ======================= COMPARISON RNEA/MAP =========================
        if(plotMAPtorques)    
            %% Comparing RNEA/MAP torques

            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;
            title(sprintf('Subject: %d, Trial: %d',subjectID, trialID))

            %RNEA
            plot1 = plot(dataTime,resRNEA.tau(:,1), 'lineWidth',2.5); hold on;
            set(plot1,'color',[1 0 0]);
            plot2 = plot(dataTime,resRNEA.tau(:,2), 'lineWidth',2.5); hold on;
            set(plot2,'color',[0 0.498039215803146 0]);

            %MAP
            plot3 = plot(dataTime,resMAP.tau_ankle, 'lineWidth',1.5,'LineStyle','--'); hold on;
            set(plot3,'color',[0 0 0]);
            plot4 = plot(dataTime,resMAP.tau_hip, 'lineWidth',1.5,'LineStyle','--'); hold on;
            set(plot4,'color',[0 0 0]);

            leg = legend('$\tau_{1,RNEA}$','$\tau_{2,RNEA}$','$\tau_{MAP}$','Location','southeast');
            set(leg,'Interpreter','latex');
            set(leg,'FontSize',18);
            xlabel('Time [s]','FontSize',20);
            ylabel('Torque [Nm]','FontSize',20);
            axis tight;
            grid on;

            if(plotPaper)
                save2pdf(fullfile(figFolder, 'MAPcomparison.pdf'),fig,600);
            end
        end
    
        %% Plot covariances
        
        sigma.Sd = myMAP.Sd;
        sigma.Sy_inv = full(mySens.sensorsParams.Sy_inv);
        sigma.SD_inv = full(dmodel.Sv_inv.matrix);
        sigma.Sprior_inv = full(dmodel.Sw_inv.matrix);
        
        if(plotCovariances)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;

            subplot(221);
            imagesc(sigma.Sd);
            colorbar;
            axis tight;
            title('Sigma (d|y)');

            subplot(222);
            imagesc(sigma.Sy_inv);
            colorbar;                                                
            axis tight;
            title('Sigma (y)');

            subplot(223);
            imagesc(sigma.SD_inv);
            colorbar;
            axis tight;
            title('Sigma (D)');

            subplot(224);
            imagesc(sigma.Sprior_inv);
            colorbar;
            axis tight;
            title('Sigma (d)');

            axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
            text(0.5, 0.99,(sprintf('Covariances (Subject: %d, Trial: %d)',subjectID, trialID)),'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',14);
        end

        %% Organising into a structure    
        finalResults(subjectID,trialID).resMAP = resMAP;
        
        clear res;

    end
    fprintf('\n');
end
%% storing results
save('./experiments/humanFixedBase/intermediateDataFiles/finalResults.mat','finalResults');

fprintf('---------\n');
fprintf('Done!\n');

