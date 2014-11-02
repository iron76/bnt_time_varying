load sim200ms.mat

close all
hh = figure
h = shadedErrorBar(time, legRightqHat(3,:)'+pi/2, sqrt(legRightqSigma(3,:))', {'k' , 'LineWidth', 2}, 0);
xlabel('time [sec]', 'FontSize', 20)
ylabel('q_3 [rad]', 'FontSize', 20)
h = legend([h.mainLine], 'estimation', 'Location', 'Southeast')
set(h, 'FontSize', 20)
grid
set(gca,'FontSize',20)
curr_ax = axis;
curr_ax(4) = curr_ax(4) +0.4;
axis(curr_ax)
print(hh, '-dpdf', 'q3Hat.pdf')

hh = figure
h2 = plot(time, legRight_q(:,[1])*pi/180+pi/2, 'k', 'LineWidth', 2)
xlabel('time [sec]', 'FontSize', 20)
ylabel('q_3 [rad]', 'FontSize', 20)
h = legend([h2], 'measure', 'Location', 'Southeast')
set(h, 'FontSize', 20)
grid
set(gca,'FontSize',20)
axis(curr_ax)
print(hh, '-dpdf', 'q3.pdf')

hh = figure
h2 = plot(time, cumsum(legRight_dq(:,1)*pi/180 + footInertial(:,9))*dtime, 'k', 'LineWidth', 2)
xlabel('time [sec]', 'FontSize', 20)
ylabel('q_3 [rad]', 'FontSize', 20)
h = legend([h2], 'simple estimation', 'Location', 'Southeast')
set(h, 'FontSize', 20)
grid
set(gca,'FontSize',20)
axis(curr_ax)
print(hh, '-dpdf', 'q3Simple.pdf')

hh = figure
h2 = plot(time(10:end-10), footInertial(10:end-10,9), 'k', 'LineWidth', 2)
xlabel('time [sec]', 'FontSize', 20)
ylabel('\omega_{4,z} [rad]', 'FontSize', 20)
h = legend([h2], 'measure', 'Location', 'Southeast')
set(h, 'FontSize', 20)
grid
set(gca,'FontSize',20)
print(hh, '-dpdf', 'omega4.pdf')


hh = figure
h = shadedErrorBar(time, legRightdqHat(3,:)', sqrt(legRightdqSigma(3,:))', {'k' , 'LineWidth', 2}, 0);
xlabel('time [sec]', 'FontSize', 20)
ylabel('dq_3 [rad/sec]', 'FontSize', 20)
h = legend([h.mainLine], 'estimation', 'Location', 'Southeast')
set(h, 'FontSize', 20)
grid
set(gca,'FontSize',20)
print(hh, '-dpdf', 'dq3Hat.pdf')
curr_ax = axis;

hh = figure
h2 = plot(time, legRight_dq(:,[1])*pi/180, 'k', 'LineWidth', 2)
xlabel('time [sec]', 'FontSize', 20)
ylabel('dq_3 [rad/sec]', 'FontSize', 20)
h = legend([h2], 'measure', 'Location', 'Southeast')
set(h, 'FontSize', 20)
grid
set(gca,'FontSize',20)
axis(curr_ax)
print(hh, '-dpdf', 'dq3.pdf')






