
clc; clear; close all;

%% selected subjects and trials
subjectList = 1:12;
trialList = 1:4;  


for subjectID = subjectList
    fprintf('\n---------\nSubject : %d\nTrial : ',subjectID);
    for trialID = trialList
         fprintf('%d, ',trialID);
        
    %% load BERDY compatible data and sensor link transforms
    load('./experiments/humanFixedBase/intermediateDataFiles/BERDYFormattedSensorData.mat');
    load('./experiments/humanFixedBase/intermediateDataFiles/sensorLinkTransforms.mat');
   
    currentTrial = BERDYFormattedSensorData(subjectID,trialID); 
    data = currentTrial.data;
    dataTime = currentTrial.data.dataTime;
    q = currentTrial.data.q';
    dq = currentTrial.data.dq';
    ddq = currentTrial.data.ddq';
    len = length(dataTime);

    currentTrialSens = sensorLinkTransforms(subjectID,trialID);
    X_imu_2 = currentTrialSens.X_imu_2;
    XStar_fp_0 = currentTrialSens.XStar_fp_0;
    XStar_0_1 = currentTrialSens.XStar_0_1;
    
    
    %%=====structure from files
    data.parts    = {'leg'         ,'torso'};
    data.labels   = {'fts'         ,'imu'  };
    data.ndof     = {6             ,6      };
    data.index    = {'1:6'         ,'1:6'  };

    %%=====structure of sensors
    sens.parts    = {'leg'         ,'torso'};       %force of the forceplate is trasmitted into the leg
    sens.labels   = {'fts'         ,'imu'  };
    sens.ndof     = {6             ,6      };

    label_to_plot = {'fts'         ,'imu'  };
  
    %% build models
    load('./experiments/humanFixedBase/intermediateDataFiles/humanThreeLinkModelFromURDF.mat');
    
    currentModel = humanThreeLinkModelFromURDF(subjectID);
   
    dmodel = currentModel.dmodel;                       %deterministic model 
    ymodel  = humanThreeLinkSens(dmodel, sens);  
   
    dmodel  = autoTreeStochastic(dmodel, 1e-5, 1e-4);   % probabilistic model for D equation (added Sv and Sw)
    ymodel  = humanThreeLinkSensStochastic(ymodel);     % probabilistic model for Y(q,dq) d = y (added Sy)
   
    myModel = model(dmodel);
    mySens  = sensors(ymodel);  
   
    myMAP  = MAP(myModel, mySens);
   
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

    Ymatrix = zeros (ymodel.m,26*dmodel.NB); 

    % The next line is time varying so must be moved inside loop.
    %Ymatrix(1:6,13:18) = sensorLinkTransforms.XStar_fp_0 * sensorLinkTransforms.XStar_0_1;

    Ymatrix(10:12,27:32) = X_imu_2(4:6,:);
    Ymatrix(13:18,20:25) = eye(6);
    Ymatrix(19:24,46:51) = eye(6);
    Ymatrix(25,26) = eye(1);
    Ymatrix(26,52) = eye(1);

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

    a_grav = [0;0;0;0;0;-9.8100]; %Featherstone-like notation

    I_c = [footIxx     0        0   ; 
              0     footIyy     0   ; 
              0        0     footIzz]; 
          
    I_0 = createSpatialInertia(I_c,footMass,posP_0);

    b_Y(1:6,1:len)   = repmat((-XStar_fp_0 * I_0 * a_grav),1,len);
    
    for i = 1 : len   
        A =R_imu_2*v{i,2}(1:3,1);
        B =((X_imu_2(4:6,1:3)*v{i,2}(1:3,1))+(R_imu_2*v{i,2}(4:6,1)));
        b_Y(10:12,i) = cross(A,B);
    end 

    clear A;
    clear B;
    %% Computing MAP method

    for i = 1 : len
    
        % reseting the only row in Ymatrix that is time varying
        Ymatrix(1:6,13:18) = XStar_fp_0 * XStar_0_1{i};
    
        myMAP = myMAP.setState(data.q(i,:)', data.dq(i,:)');
        myMAP = myMAP.setY(data.y(:,i));
        myMAP = myMAP.setYmatrix(Ymatrix);
        myMAP = myMAP.setBias(b_Y(:,i));
        myMAP = myMAP.solveID();
  
        res.d(:,i)    = myMAP.d;
        res.Sd(:,:,i) = myMAP.Sd; %full() passing from sparse to double matrix
        res.y(:,i)    = (Ymatrix * res.d(:,i)) + b_Y(:,i); 
        res.Sy(:,:,i) = Ymatrix * res.Sd(:,:,i) * Ymatrix';
    
         if mod(i-1,100) == 0
                fprintf('Processing %d %% of the dataset\n', round(i/len*100));
         end
     end
% 
% 
% %====plot of Ymatrix
% imagesc(Ymatrix)
% colorbar
% title('Y matrix','FontSize',15);

% ========end MAP
%% Plot overlapped plots
% 
% py = [0; cumsum(cell2mat(myMAP.IDsens.sensorsParams.sizes))];
% for l = 1 : length(label_to_plot)
%    for k = 1 : myMAP.IDsens.sensorsParams.ny
%       if strcmp(myMAP.IDsens.sensorsParams.labels{k}, label_to_plot{l})
%          figure
%          J = myMAP.IDsens.sensorsParams.sizes{k};
%          I = py(k)+1 : py(k)+J;
%          colors = ['r', 'g', 'b'];
%          for j = 1 : J
%             subplot(2, ceil(J/2), j)
%              hold on;
%             shadedErrorBar(data.time, data.y(I(j),:), sqrt(data.Sy(I(j), :)), {[colors(mod(j,3)+1) '--'] , 'LineWidth', 1}, 0);
%             plot(data.time, data.y(I(j),:), colors(mod(j,3)+1) , 'LineWidth', 1);
% 
%             title(strrep(myMAP.IDsens.sensorsParams.labels{k}, '_', '~'));
%          end
%       end
%    end
% end
%%
% % Plotting covariance matrices (some examples)
% 
% %===Covariance Sigma(d|y)
% figure
% subplot(2,1,1)
% imagesc(res.Sd(:,:,1))
% colorbar
% title('Covariance Sigma(d|y) matrix','FontSize',15);
% 
% %===Covariance Sigma(y)
% subplot(2,1,2)
% imagesc(res.Sy(:,:,1))
% colorbar
% title('Covariance Sigma(y) matrix','FontSize',15);


    %% Rerrange solution
 
    % d  = zeros(26*dmodel.NB, len);
    % Sd = zeros(26*dmodel.NB, 26*dmodel.NB, len);

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
    %% organising MAP results
       MAPresults(subjectID,trialID).res = res;

    %% ======================= LEAST SQUARE SOLUTION ==========================
    %  
    % %Using RNEA class for fetting D and b_D 
    % ymodel_RNEA  = autoSensRNEA(dmodel);
    % mySens_RNEA  = sensors(ymodel_RNEA);
    % myRNEA       = RNEA(myModel, mySens_RNEA);
    %       
    % %Ordering y_RNEA in the form [fx1 fx2 ddq1 ddq2]
    % y_RNEA_f = zeros(6*dmodel.NB, len);
    % y_RNEA_ddq = zeros(dmodel.NB, len);
    % fx = cell(dmodel.NB);
    %   
    % for i = 1 : dmodel.NB
    %    for t = 1 : len
    %           fx{i,1} = zeros(6,1); 
    %           y_RNEA_f(((1:6)+(i-1)*6), t) = [fx{i,1}];
    %           y_RNEA_ddq(i, t) = [data.ddq(i,t)];
    %    end
    %    y_RNEA = [y_RNEA_f ; y_RNEA_ddq];
    % end
    % 
    % %computing d_RNEA
    % d_RNEA = zeros (26*myRNEA.IDmodel.modelParams.NB,len);
    %    for i = 1 : len
    %         myRNEA = myRNEA.setState(data.q(:,i), data.dq(:,i));
    %         myRNEA = myRNEA.setY(y_RNEA(:,i));
    %         myRNEA = myRNEA.solveID();
    %        
    %         d_RNEA(:,i) = myRNEA.d; 
    %         b_RNEA(:,i) = myRNEA.b.matrix;
    %    end
    %    
    % % since:       | D |     | b_D |   | 0 |
    % %              |   | d + |     | = |   |
    % %              | Y |     | b_Y |   | y |
    % %                       
    % %                   D_invDyn   b_invDyn  y_invDyn
    % 
    % 
    %    D_invDyn = [myRNEA.D.matrix;
    %                    Ymatrix    ];
    %   
    %    y_d  = zeros(38,len);
    %    d_ls =  zeros(size(d_RNEA));
    %           
    %    for k = 1:len
    %    
    %        b_invDyn(:,k) = [b_RNEA(:,k);
    %                           b_Y(:,k) ];
    %        
    %        y_invDyn(:,k) = [  y_d(:,k)  ;
    %                         data.y(:,k)];
    %      
    %             
    %        %least square solution of D_invDyn*d - RS = 0
    %        RS =  y_invDyn - b_invDyn;  
    %        d_ls(:,k) = pinv(D_invDyn)* RS(:,k); %equivalent to  d_ls_temp(:,k)= D_invDyn\RS(:,k);  
    %                     
    %    end
    % 
    % % ========end LS
    
    %% ======================= COMPARISON RNEA/MAP/LS =========================
    %%
    
    q = currentTrial.data.q;
    dq = currentTrial.data.dq;
    ddq = currentTrial.data.ddq;
    
    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;
    
    subplot(311);
    plot1 = plot(dataTime,(180/pi)*q(:,1),'lineWidth',2.0); hold on;
    set(plot1,'color',[1 0 0]);
    plot2= plot(dataTime,(180/pi)*q(:,2),'lineWidth',2.0); hold on;
    set(plot2,'color',[0 0.498039215803146 0]);
    leg = legend('$q_1$','$q_2$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Joint Angle [deg]','FontSize',18);
    axis tight;
    grid on;

    subplot(312);
    plot1 = plot(dataTime,(180/pi)*dq(:,1),'lineWidth',2.0); hold on;
    set(plot1,'color',[1 0 0]);
    plot2= plot(dataTime,(180/pi)*dq(:,2),'lineWidth',2.0); hold on;
    set(plot2,'color',[0 0.498039215803146 0]);
    leg = legend('$\dot q_1$','$\dot q_2$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Joint Velocity [deg/s]','FontSize',18);
    axis tight;
    grid on;

    subplot(313);
    plot1 = plot(dataTime,tau(:,1), 'lineWidth',2.0); hold on;
    set(plot1,'color',[1 0 0]);
    plot2 = plot(dataTime,tau(:,2), 'lineWidth',2.0); hold on;
    set(plot2,'color',[0 0.498039215803146 0]);
    leg = legend('$\tau_1$','$\tau_2$','Location','northeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',18);
    ylabel('Torque [Nm]','FontSize',18);
    axis tight;
    grid on;
    
    %% Comparing RNEA/MAP torques

    fig = figure();
    axes1 = axes('Parent',fig,'FontSize',16);
    box(axes1,'on');
    hold(axes1,'on');
    grid on;

    % RNEA
    plot1 = plot(dataTime,tau(:,1), 'lineWidth',2.5); hold on;
    set(plot1,'color',[1 0 0]);
    plot2 = plot(dataTime,tau(:,2), 'lineWidth',2.5); hold on;
    set(plot2,'color',[0 0.498039215803146 0]);

    % MAP
    plot3 = plot(dataTime,res.tau_ankle, 'lineWidth',1.5,'LineStyle','--'); hold on;
    set(plot3,'color',[1 0 0]);
    plot4 = plot(dataTime,res.tau_hip, 'lineWidth',1.5,'LineStyle','--'); hold on;
    set(plot4,'color',[0 0.498039215803146 0]);

    % % LS
    % plot5 = plot(data.time,d_ls(19,:), 'lineWidth',1.5,'LineStyle',':'); hold on;
    % set(plot5,'color',[1 0 0]);
    % plot6 = plot(data.time,d_ls(45,:), 'lineWidth',1.5,'LineStyle',':'); hold on;
    % set(plot6,'color',[0 0.498039215803146 0]);

    %leg = legend('$\tau_{1,RNEA}$','$\tau_{2,RNEA}$','$\tau_{2,MAP}$','$\tau_{2,MAP}$','$\tau_{1,LS}$','$\tau_{2,LS}$','Location','southeast');
    leg = legend('$\tau_{1,RNEA}$','$\tau_{2,RNEA}$','$\tau_{2,MAP}$','$\tau_{2,MAP}$','Location','southeast');
    set(leg,'Interpreter','latex');
    set(leg,'FontSize',18);
    xlabel('Time [s]','FontSize',20);
    ylabel('Torque [Nm]','FontSize',20);
    axis tight;
    grid on;

    end
    fprintf('\n');
end
%% storing results
save('./experiments/humanFixedBase/intermediateDataFiles/MAPresults.mat','MAPresults');



% %% Simulate output in MAP
% 
%  for  ind = 1:26
%      
%         figure;
%                 
%        
%         %====== Comparison between MAP prediction and actual data
%         
%         %== simulate output in MAP
%         y_pred_MAP= myMAP.simY(res.d);     % without b_Y
%         y_pred_MAP = y_pred_MAP + b_Y;     % adding b_Y
%         
%         subplot(2,1,1);
%         plot1 = plot(data.time,y_pred_MAP(ind,:), 'lineWidth',1.0, 'LineStyle','--'); hold on;
%         set(plot1,'color',[1 0 0]);
%         plot2 = plot(data.time,data.y(ind,:), 'lineWidth',1.0); hold on;
%         set(plot2,'color',[0 0 1]);
% 
%         leg = legend('MAP Pred','Input data','Location','northeast');
%         %set(leg,'Interpreter','latex');
%         set(leg,'FontSize',15);
%         xlabel('Time [s]','FontSize',15);
%         %ylabel('Torque[Nm]','FontSize',20);
%         
%         if (ind >=1 && ind<4) 
%         title(sprintf('Figure %d : f1-moments ',ind));
%         elseif (ind >=4 && ind<7) 
%         title(sprintf('Figure %d : f1-forces ',ind));
%         elseif (ind >=7 && ind<10)
%         title(sprintf('Figure %d : a2-ang',ind));
%         elseif (ind >=10 && ind<13)
%         title(sprintf('Figure %d : a2-lin ',ind));
%         elseif (ind >=13 && ind<19)
%         title(sprintf('Figure %d : fx1 ',ind));
%         elseif (ind >=19 && ind<25)
%         title(sprintf('Figure %d : fx2 ',ind));
%         elseif (ind ==25)
%         title(sprintf('Figure %d : ddq1 ',ind));
%         elseif (ind ==26)
%         title(sprintf('Figure %d : ddq2 ',ind));
%         end
%         
%         axis tight;
%         grid on;
%         
%         
%         %====== Comparison between RNEA prediction and actual data
%         
%         %== simulate output in RNEA
%         y_pred_RNEA = myMAP.simY(d);
%          
%         subplot(2,1,2);
%         plot1 = plot(data.time,y_pred_RNEA(ind,:), 'lineWidth',1.0, 'LineStyle','--'); hold on;
%         set(plot1,'color',[1 0 0]);
%         plot2 = plot(data.time,data.y(ind,:), 'lineWidth',1.0); hold on;
%         set(plot2,'color',[0 0 1]);
%         
%         leg = legend('RNEA Pred', 'Input data','Location','northeast');
%         %set(leg,'Interpreter','latex');
%         set(leg,'FontSize',15);
%         xlabel('Time [s]','FontSize',15);
%         %ylabel('Torque[Nm]','FontSize',20);
%         %title(sprintf('Figure %d',ind));
%         axis tight;
%         grid on;
% 
%  end
% 
%  %% ============================= tests ==================================
%  load ('sensorLinkTransforms.mat');
%  %% test 1 : comparison between a_2 measured by sensor and a_2 in vector d_RNEA
%  
%  %  
%  % -in Featherstone notation:              | a_2,ang |     |      domega_2      |
%  %                                   a_2 = |         |   = |                    |   
%  %                                         | a_2,lin |     | ddr_2 - omega x dr |
%  %
%  % -a_2 measured is in imu frame;
%  % -a_2 in vector d_RNEA is in frame associated to link2;
%  % -both a_2 are compared in frame associated to link2
% 
%  
%  a_imu_2real = zeros (6,len);
%  a_2_2real   = zeros (6,len);
% 
%  for i = 1:len
%     % in sensor frame
%     a_imu_2real(:,i) = data.y(7:12,i) - b_Y(7:12,i); 
%  
%     % in frame associated to link2 
%     a_2_2real(:,i) = (X_imu_2)'* a_imu_2real(:,i);
%  end
%  
%  % now we can compare a_2_2real measured of sensor (transformed in link
%  % frame) with a_2_2 of d vector.
%  
%   %% test 2 : comparison between f measured by sensor and f in vector d_RNEA  
%  
%  % -f_1 measured is in fp frame;
%  % -f_1 in vector d_RNEA is in frame associated to link1;
%  % -both f_1 are compared in frame associated to link1 
%  
% f_fp_1real = length (a_imu_2real);
% f_1_1real = length (a_imu_2real);
% 
%  % computing 1_XStar_fp :   1_XStar_fp = 1_XStar_0 * 0_XStar_fp
%   XStar_1_fp = (XStar_0_1)' * (XStar_fp_0)';
% 
% 
% for i = 1:len
%     % in sensor frame
%     f_fp_1real(:,i) = data.y(1:6,i) + b_Y(1:6,i);
%     
%     % in frame associated to link1 
%     f_1_1real(:,i) = XStar_1_fp * f_fp_1real(:,i);
% end
%  
%  % now we can compare f_1_1real measured of sensor (transformed in link
%  % frame) with f_1_1 of d vector.
%  
% 
