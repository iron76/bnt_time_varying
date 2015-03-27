clear all
close all
clc
load preprocess2.mat
load d0.mat

[~, n]  = size(y);
NB      = myModel.modelParams.NB;
m       = 2;
dx      = 0;

%        'torso_pitch' 'torso_roll' 'torso_yaw'
Q  = [              1,           1,          1];
Sx = [              1,           1,          1];

%         'l_hip_pitch' 'l_hip_roll' 'l_hip_yaw' 'l_knee' 'l_ankle_pitch' 
Q  = [Q,             1,           1,          1,       0,              1];
Sx = [Sx,            1,           1,          1,    1e-5,              1];

%         'l_shoulder_pitch' 'l_shoulder_roll' 'l_shoulder_yaw' 'l_elbow' 'l_wrist_prosup' 'l_ankle_roll'
Q  = [ Q,                 1,                1,               1,        1,               1,              1];
Sx = [Sx,                 1,                1,               1,        1,               1,              1];

%         'r_hip_pitch' 'r_hip_roll' 'r_hip_yaw' 'r_knee' 'r_ankle_pitch'
Q  = [Q,             1,           1,          1,       0,              1];
Sx = [Sx,            1,           1,          1,    1e-5,              1];

%         'r_shoulder_pitch' 'r_shoulder_roll' 'r_shoulder_yaw' 'r_elbow' 'r_wrist_prosup' 'r_ankle_roll'
Q  = [Q,                  1,                1,               1,        1,               1,              1];
Sx = [Sx,                 1,                1,               1,        1,               1,              1];

Q        = 1e-4*diag([Q, Q]);
Sq       = 0.4   * pi/180;
Sdq      = 0.1   * pi/180;
Sx0      = 1e2*diag([Sx*Sq, Sx*Sdq]);

%% iCub model
ymdl    = iCubSens(dmodel, sens);
ymodel  = iCubSensDNEA(dmodel, ymdl, sens, mask_q, mask_dq);
dmodel  = autoTreeStochastic(dmodel, 1e-2, 1e4);
ymodel  = iCubSensStochastic(ymodel);
myModel = model(dmodel);
mySens  = sensors(ymodel);

%% Compute solution
myDNEA    = DNEA(myModel, mySens);

res.d  = zeros(26*NB, n);
res.Sd = zeros(26*NB, 26*NB, n);
res.x  = zeros(2 *NB, n);
res.Sx = zeros(2 *NB, 2* NB, n);

for j = 1 : m
   for i = 1 : n
      x      = [data.q(:,i); data.dq(:,i)];
      x_pri  = x + dx;
      
      myDNEA = myDNEA.setState(x_pri(1:NB ,1),x_pri(NB+1:end,1));
      myDNEA = myDNEA.setY(data.y(:,i));
      
      myDNEA = myDNEA.setD(d0);
      myDNEA = myDNEA.setDprior(d0);
      myDNEA = myDNEA.setXprior(x_pri);
      myDNEA = myDNEA.setXvariance(Sx0);
      myDNEA = myDNEA.solveID();
      
      Y = cell2mat(myDNEA.IDsens.sensorsParams.Y);
      
      res.d(:,i)    = myDNEA.d;
      res.x(:,i)    = myDNEA.x;
      res.Sd(:,:,i) = myDNEA.Sd;
      res.Sx(:,:,i) = myDNEA.Sx;
      
      res.y(:,i)    = Y * [res.d(:,i); res.x(:,i)];
      res.Sy(:,:,i) = Y * blkdiag(res.Sd(:,:,i), res.Sx(:,:,i))* Y';
      if mod(i-1,10) == 0
         fprintf('Processing %d %% of the dataset\n', round(i/n*10)*10);
      end
      d0  = res.d(:,i);
      dx  = dx + myDNEA.x - x_pri;
      Sx0 = myDNEA.Sx + Q;
   end
   data.dx_hat(:,j) = dx;
   
   Sx0 = 1e2*diag([Sx*Sq, Sx*Sdq]);
   dx  = Sx0 * randn(2*NB,1)* pi/180;
   
   
   %% Rerrange solution
   for i = 1 : NB
      for j = 1 : n
         link = strrep(myDNEA.IDmodel.modelParams.linkname{i}, '+', '_');
         di   = ['res.d_'   link  '(:,j)'];
         ind  = '1 + 26*(i-1) : 26*(i-1) + 26';
         eval([di '   = res.d(' ind ',j);'])
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
         eval(['res.tau_'  link '(:,j)   =  res.d(' ind '        ,j);'])
         eval(['res.Stau_' link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
         %fx
         ind  = '20 + 26*(i-1) : 26*(i-1) + 25';
         eval(['res.fx_'   link '(:,j)   =  res.d(' ind '        ,j);'])
         eval(['res.Sfx_'  link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
         %d2q
         ind  = '26 + 26*(i-1) : 26*(i-1) + 26';
         eval(['res.d2q_'  link '(:,j)   =  res.d(' ind '        ,j);'])
         eval(['res.Sd2q_' link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
      end
   end
   
   res.q  = res.x(    1:NB, :);
   res.dq = res.x(NB+1:end, :);
   
   %% Plot overlapped plots
   % py = [0; cumsum(cell2mat(myDNEA.IDsens.sensorsParams.sizes))];
   % label_to_plot = {'imu'                                      , 'la_acc'                 , 'lf_acc'                    , 'lh_acc'                              , 'ra_acc'                  , 'rf_acc'                    , 'rh_acc'                              , 'to_acc'                                  , 'la_fts'                    , 'ra_fts'                   , 'll_fts'               , 'rl_fts'               , 'lf_fts'                    , 'rf_fts'                     , 'lh_gyr'                               , 'rh_gyr'                              };
   % for l = 1 : length(label_to_plot)
   %    for k = 1 : myDNEA.IDsens.sensorsParams.ny - dmodel.NB
   %       if strcmp(myDNEA.IDsens.sensorsParams.labels{k}, label_to_plot{l})
   %          figure
   %          J = myDNEA.IDsens.sensorsParams.sizes{k};
   %          I = py(k)+1 : py(k)+J;
   %          colors = ['r', 'g', 'b'];
   %          for j = 1 : J
   %             subplot(2, ceil(J/2), j)
   %             hold on;
   %             shadedErrorBar(data.time, data.y(I(j),1:length(data.time)), sqrt(data.Sy(I(j), 1:length(data.time))), {[colors(mod(j,3)+1) '--'] , 'LineWidth', 1}, 0);
   %             shadedErrorBar(data.time,  res.y(I(j),1:length(data.time)), sqrt(reshape(res.Sy(I(j), I(j), 1:length(data.time)), 1, n)), {colors(mod(j,3)+1) , 'LineWidth', 2}, 0);
   %
   %             title(strrep(myDNEA.IDsens.sensorsParams.labels{k}, '_', '~'));
   %          end
   %       end
   %    end
   % end
   
   rad2deg = 180/pi;
   for i = 1 : NB
      figure(i)
      subplot(211)
      hold on
      h1 = shadedErrorBar(data.time, rad2deg*data.q(i,1:length(data.time))', sqrt(ones(1,n)*Sx0(i)), {'k--' , 'LineWidth', 1}, 0);
      h2 = shadedErrorBar(data.time,  rad2deg*res.q(i,1:length(data.time))', sqrt(reshape(res.Sx(i, i, 1:length(data.time)), 1, n)), {'k' , 'LineWidth', 2}, 1);
      grid
      ax = gca;
      set(ax, 'FontSize', 24)
      legend('q', 'q_{map}')
      ylabel('q [deg]')
      [legh,objh,outh,outm] =legend([h1.mainLine h2.mainLine], 'q', 'q_{map}', 'Location', 'Southeast');
      title(strrep(myDNEA.IDmodel.modelParams.jointname{i}, '_', '-'));
      set(legh,'linewidth',4);
      
      subplot(212)
      hold on
      plot(data.time,  rad2deg*(res.q(i,1:length(data.time))'- data.q(i,1:length(data.time))'), 'k', 'LineWidth', 2)
      grid
      ax = gca;
      set(ax, 'FontSize', 24)
      xlabel('time[s]')
      ylabel('e [deg]')
      eval(['print -dpdf ' strrep(myDNEA.IDmodel.modelParams.jointname{i}, '_', '-') '.pdf'])
      
      fprintf('Joint %s mean: %f, variability: %f \n', myDNEA.IDmodel.modelParams.jointname{i}, rad2deg*mean(data.dx_hat(i,:)), rad2deg*std(data.dx_hat(i,:)))
   end
end


