% load data
clear;clc;

load('preProcessedSensorData.mat');
load('humanThreeLinkModelFromURDF.mat');


q1 = processedSensorData.q1;
q2 = processedSensorData.q2;

dq1 = zeros (length(q1),1);
dq2 = zeros (length(q2),1);

ddq1 = zeros (length(q1),1);
ddq2 = zeros (length(q2),1);


q  = [q1,q2];
dq  = [dq1,dq2];
ddq  = [ddq1,ddq2];


%tau = zeros(size(q));
%a = cell (size(q));
for i = 1:size(q)
      [tau_i, a_i, fB_i, f_i] = ID( humanThreeLink_dmodel, q(i,:), dq(i,:), ddq(i,:));
      tau(i,:) = tau_i';
%      a(i,:) = a_i';
%     fB(i,:) = fB_i';
%     f(i,:) = f_i';
end


figure;
plot(tau); hold on;
xlabel('time (sec)');
ylabel('torque (Nm)');
legend('\tau_1','\tau_2');
axis tight;
