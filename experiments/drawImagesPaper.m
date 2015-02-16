clear all
close all
load predict.mat

h1 = figure;
XYZ = {'x', 'y', 'z'};
for i = 1 : 3
   subplot(3,1,i);
   hf = shadedErrorBar(data.time,  res.fx_r_upper_foot_r_foot_r_sole(i,:), res.Sfx_r_upper_foot_r_foot_r_sole(i,:), {'k--' , 'LineWidth', 2}, 0);
   xlabel('time [sec]', 'FontSize', 20);
   ylabel(['f_' XYZ{i} '[N]'], 'FontSize', 20);
   set(gca,'FontSize',20);
   % plot(data.time,  res.fx_r_upper_foot_r_foot_r_sole(i,:), 'k--', 'LineWidth', 2);
   grid on
end

h2 = figure;
for i = 4 : 6
   subplot(3,1,i-3);
   hu = shadedErrorBar(data.time,  res.fx_r_upper_foot_r_foot_r_sole(i,:), res.Sfx_r_upper_foot_r_foot_r_sole(i,:), {'k--' , 'LineWidth', 2}, 0);
   xlabel('time [sec]', 'FontSize', 20);
   ylabel(['\mu_' XYZ{i-3} '[N]'], 'FontSize', 20);
   set(gca,'FontSize',20)
   % plot(data.time,  res.fx_r_upper_foot_r_foot_r_sole(i,:), 'k--', 'LineWidth', 2);
   grid on
end

load predictLearnt.mat

figure(h1)
for i = 1 : 3
   subplot(3,1,i);
   hold on
   colors = ['r', 'g', 'b'];   
   % shadedErrorBar(data.time,  res.fx_r_upper_foot_r_foot_r_sole(i,:), res.Sfx_r_upper_foot_r_foot_r_sole(i,:), {colors(mod(i,3)+1) , 'LineWidth', 2}, 0);
   h = plot(data.time,  res.fx_r_upper_foot_r_foot_r_sole(i,:), 'k', 'LineWidth', 2);
   grid on
end
h = legend([h hf.mainLine], 'After learning', 'Before learning', 'Location', 'Southeast');


figure(h2)
for i = 4 : 6
   subplot(3,1,i-3);
   hold on
   colors = ['r', 'g', 'b'];   
   % shadedErrorBar(data.time,  res.fx_r_upper_foot_r_foot_r_sole(i,:), res.Sfx_r_upper_foot_r_foot_r_sole(i,:), {colors(mod(i,3)+1) , 'LineWidth', 2}, 0);
   h = plot(data.time,  res.fx_r_upper_foot_r_foot_r_sole(i,:), 'k', 'LineWidth', 2);
   grid on
end
h = legend([h hu.mainLine], 'After learning', 'Before learning', 'Location', 'Southeast')


print(h1, '-dpdf', 'forces.pdf')
print(h2, '-dpdf', 'torques.pdf')

clear all
close all
load learn.mat

h3 = figure;
subplot(211)
h = plot(real(ll), 'k', 'LineWidth', 2);
xlabel('iterations', 'FontSize', 20)
ylabel(['Log-likelihood'], 'FontSize', 20)
set(gca,'FontSize',20)
grid
print(h3, '-dpdf', 'learn.pdf')


