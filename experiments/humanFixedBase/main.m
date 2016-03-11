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
        % Computing d using Newton-Euler with inverse dynamics.
        % Function ID was modified for having velocities.  So the path is not:
        % ../bnt_time_varying/extern/featherstone/dynamics/ID.m
        % but
        % ../bnt_time_varying/experiments/humanFixedBase/helperFunction/IDv.m

        tau = zeros(size(q))';
        a = cell (size(q))';
        fB = cell (size(q))';
        f = cell (size(q))';
        fx = zeros (6,1);

        v = cell(size(q))';

        d_temp = zeros(26*dmodel.NB,1);
        d = zeros (26*dmodel.NB, len);

        fext    = cell(1,2);
        for i = 1 : dmodel.NB
             fext{i}    = fx;
        end

        for i = 1:len

            [tau_i, a_i,v_i,fB_i, f_i] = IDv(dmodel, q(:,i), dq(:,i), ddq(:,i), fext);
            tau(i,:) = tau_i;
            a(i,:) = a_i;
            fB(i,:) = fB_i;
            f(i,:) = f_i;  
            v(i,:) = v_i;

            for j = 1 : dmodel.NB
                 d_temp((1:26)+(j-1)*26) = [a_i{j}; fB_i{j}; f_i{j}; tau(i,j); fx; ddq(j,i)];
            end

            d(:,i) = d_temp;
        end

        clear d_temp;
        clear tau_i; 
        clear a_i;
        clear fB_i;
        clear f_i;
        clear v_i;
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

        for i = 1 : len   
            A =R_imu_2*v{i,2}(1:3,1);
            B =((X_imu_2(4:6,1:3)*v{i,2}(1:3,1))+(R_imu_2*v{i,2}(4:6,1)));
            b_Y(10:12,i) = cross(A,B);
        end 

        clear A;
        clear B;
        %% Computing MAP method

        for i = 1 : len

            myMAP = myMAP.setState(data.q(i,:)', data.dq(i,:)');
            myMAP = myMAP.setY(data.y(:,i));
            myMAP = myMAP.setYmatrix(Ymatrix{i});
            myMAP = myMAP.setBias(b_Y(:,i));
            myMAP = myMAP.solveID();

            res.d(:,i)       = myMAP.d;
            res.Sd(:,:,i)    = myMAP.Sd; %full() passing from sparse to double matrix
            res.Ymatrix{i,1} = myMAP.IDsens.sensorsParams.Y; 
            res.b_Y(:,i)     = myMAP.IDsens.sensorsParams.bias;
            %res.y(:,i)      = (Ymatrix{i} * res.d(:,i)) + b_Y(:,i); 
            %res.Sy(:,:,i)   = Ymatrix{i} * res.Sd(:,:,i) * Ymatrix{i}';

             if mod(i-1,100) == 0
                    fprintf('Processing %d %% of the dataset\n', round(i/len*100));
             end
        end
        % 
        % %====plot of Ymatrix
        % imagesc(Ymatrix)
        % colorbar
        % title('Y matrix','FontSize',15);

        % ========end MAP
        %% Rerrange solution

        for i = 1 : dmodel.NB
             for j = 1 : len

                 link = strrep(myMAP.IDmodel.modelParams.linkname{i}, '+', '_');
                 joint = strrep(myMAP.IDmodel.modelParams.jointname{i}, '+', '_');

                 di   = ['res.d_'   link  '(:,j)'];
                 ind  = '1 + 26*(i-1) : 26*(i-1) + 26';
                 eval([di '   = d(' ind ',j);'])

                 %a
                 ind  = '1 + 26*(i-1) : 26*(i-1) +  6';
                 eval(['res.a_'    link '(:,j)   =  res.d(' ind '        ,j);'])
                 eval(['res.Sa_'   link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
                 %fB
                 ind  = '7 + 26*(i-1) : 26*(i-1) + 12';
                 eval(['res.fB_'   link '(:,j)   =  res.d(' ind '        ,j);'])
                 eval(['res.SfB_'  link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
                 %f
                 ind  = '13 + 26*(i-1) : 26*(i-1) + 18';
                 eval(['res.f_'    link '(:,j)   =  res.d(' ind '        ,j);'])
                 eval(['res.Sf_'   link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
                 %tau
                 ind  = '19 + 26*(i-1) : 26*(i-1) + 19';
                 eval(['res.tau_'  joint '(:,j)   =  res.d(' ind '        ,j);'])
                 eval(['res.Stau_' joint '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
                 %fx
                 ind  = '20 + 26*(i-1) : 26*(i-1) + 25';
                 eval(['res.fx_'   link '(:,j)   =  res.d(' ind '        ,j);'])
                 eval(['res.Sfx_'  link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
                 %ddq
                 ind  = '26 + 26*(i-1) : 26*(i-1) + 26';
                 eval(['res.d2q_'  joint '(:,j)   =  res.d(' ind '        ,j);'])
                 eval(['res.Sd2q_' joint '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
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
        plot1 = plot(dataTime,tau(:,1), 'lineWidth',2.5); hold on;
        set(plot1,'color',[1 0 0]);
        plot2 = plot(dataTime,tau(:,2), 'lineWidth',2.5); hold on;
        set(plot2,'color',[0 0.498039215803146 0]);

        %MAP
        plot3 = plot(dataTime,res.tau_ankle, 'lineWidth',1.5,'LineStyle','--'); hold on;
        set(plot3,'color',[1 0 0]);
        plot4 = plot(dataTime,res.tau_hip, 'lineWidth',1.5,'LineStyle','--'); hold on;
        set(plot4,'color',[0 0.498039215803146 0]);

        leg = legend('$\tau_{1,RNEA}$','$\tau_{2,RNEA}$','$\tau_{2,MAP}$','$\tau_{2,MAP}$','Location','southeast');
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
        
        if(plotCovariances)
            fig = figure();
            axes1 = axes('Parent',fig,'FontSize',16);
            box(axes1,'on');
            hold(axes1,'on');
            grid on;

            subplot(221);
            imagesc(full(myMAP.Sd));
            colorbar;
            axis tight;
            title('Sigma (d|y)');

            subplot(222);
            imagesc(full(mySens.sensorsParams.Sy_inv));
            colorbar;                                                
            axis tight;
            title('Sigma (y)');

            subplot(223);
            imagesc(full(dmodel.Sv_inv.matrix));
            colorbar;
            axis tight;
            title('Sigma (D)');

            subplot(224);
            imagesc(full(dmodel.Sw_inv.matrix));
            colorbar;
            axis tight;
            title('Sigma (d)');

            axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping' ,'off');
            text(0.5, 0.99,(sprintf('Covariances (Subject: %d, Trial: %d)',subjectID, trialID)),'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',14);
        end

        %% Organising into a structure    
        MAPresults(subjectID,trialID).MAPres = res;
        MAPresults(subjectID,trialID).data.y = data.y;
        MAPresults(subjectID,trialID).data.Sy = data.Sy;
        MAPresults(subjectID,trialID).Sv_inv = dmodel.Sv_inv.matrix;
        MAPresults(subjectID,trialID).Sw_inv = dmodel.Sw_inv.matrix;
        clear res;

    end
    fprintf('\n');
end
%% storing results
save('./experiments/humanFixedBase/intermediateDataFiles/MAPresults.mat','MAPresults');

fprintf('---------\n');
fprintf('Done!\n');

