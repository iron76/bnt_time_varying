
% Sestri Leveante, 23 Luglio 2014
%
% This piece of code simulates the problem of estimating dynamic quantities
% for a single rigid body with distributed force/torque measurements and
% distributed gyro/accelerometers measurements. The motion is governed by
% the following differential equation:
%
% m    dv^B    + S(omega^B) (m       v^B) = f^B_1  + ... + f^B_n + mg^B
%
% I^B domega^B + S(omega^B) (I^B omega^B) = mu^B_1 + ... + mu^B_n
%
%                                 dq      = 1/2 q_omega_B . q
%
% where we defined the following quantities:
%
% I^B    : inertia in the body reference frame
% m      : mass of the rigid body
% f^B_i  : i-th force expressed in the body reference frame
% mu^B_i : i-th torque expressed in the body reference frame
% omega^B: angular velocity expressed in the body reference frame
% v^B    : linear velocity in the body reference frame
% q      : quaternion representing the rigid body orientation

clear all
close all
clc

dt      = 0.01;     % sampling time
T       = 1.0;       % time span
sigma_f = 0.5;       % output error
sigma_u = 0.05;      % output error
sigma_a = 0.05;      % output error
n       = 16;        % state dimension
m       = 9;         % output dimension

model.I  = diag([0.05 0.02 0.03]);
model.m  = 20;
model.u  = 1;
model.v  = 0.1;
R0       = eye(3);
model.x0 = [ones(6,1); ones(6,1); dcm2q(R0)'];
model.dt = dt;
model.g  = 9.81;

% [tv1, f1]=ode45(@(t, x) rigidBodyDifferentialEquation(t, x, model),[0 5], model.x0, []);
% odeSettings = odeset('Mass', @(t,y)massMatrix(t,y,model), 'MStateDependence', 'strong');
% [tv3, f3]=ode45(@(t, x) rigidBodyDifferentialEquationImplicit(t, x, model),[0 5], model.x0, odeSettings);


[t, x] = integrateForwardDynamics(0:dt:T, model.x0, model);

% Let's compute the output of used in the Kalman filter. In this specific
% case we will use dv^B, f^B and mu^B. In practice:
%
% y = [dv^B, f^B, mu^B];
%

for i = 1:length(t)
    y(:,i) = rigidBodyOutput(t(i),x(i,:)',model);
end
y = y';

R         = diag([sigma_a.*ones(1,3), sigma_f.*ones(1,3), sigma_u.*ones(1,3)]);
y         = y + rand(size(y))*chol(R);

% Let's define the state of the Kalman filter. In this case we use an
% augmented state that includes also the applied forces and torques f^B and
% mu^B. Therefore we have:
%
% X = [dv^B, domega^b, f^B, mu^B, q] = [x, f^B, mu^B, q];
%
% Notice that since we have:
%
% dX = f(x) + Af*[f^B, mu^B]
% Y  = [f(x) + Af*[f^B, mu^B]](1:3)
%       f^B, mu^B];
%
% its discretization is:
%
% Xn     = X + df/dx(x)*dt + Af*[f^B, mu^B]*dt
% dh/dX  = [df/dx(1:3) Af(1:3)
%          0           I]

Xhat             = zeros(size(x));
eta              = zeros(size(y));
Xhat(1,:)        = rand(size(model.x0')).*100 - 50 ;
% Xhat(1, 4: 6)    = model.x0( 4: 6)';
% Xhat(1,13:16)    = model.x0(13:16)';
P(:,:,1)         = eye(n);
Q                = diag([ones(6,1).*dt*100; ones(6,1).*100; ones(4,1).*dt*.1;]);

Af        = [...
    1/model.m*eye(3), zeros(3,3);
    zeros(3,3)     inv(model.I)];

model.u  = 0;
model.v  = 0;

for i = 2:length(t)+1
    df = rigidBodyDynamicsDerivatives(...
        Xhat(i-1,1) ,Xhat(i-1,2) ,Xhat(i-1,3),...
        Xhat(i-1,4) ,Xhat(i-1,5) ,Xhat(i-1,6),...
        Xhat(i-1,7) ,Xhat(i-1,8) ,Xhat(i-1,9),...
        Xhat(i-1,10),Xhat(i-1,11),Xhat(i-1,12),...
        Xhat(i-1,13),Xhat(i-1,14),Xhat(i-1,15),Xhat(i-1,16),...
        model.I(1,1),model.I(2,2),model.I(3,3),model.m,model.g);
    A = df.*dt + eye(size(df));
    % Prediction step update
    [xn, Pn] = predictStepKF(Xhat(i-1,:)', P(:,:,i-1), A, Q, model);
    % [xn, Pn] = predictStepKF(x(i-1,:)', P(:,:,i-1), A, Q, model);
    xn = xn';
    
    dh = rigidBodyOutputsDerivatives(...
        Xhat(i-1,1) ,Xhat(i-1,2) ,Xhat(i-1,3),...
        Xhat(i-1,4) ,Xhat(i-1,5) ,Xhat(i-1,6),...
        Xhat(i-1,7) ,Xhat(i-1,8) ,Xhat(i-1,9),...
        Xhat(i-1,10),Xhat(i-1,11),Xhat(i-1,12),...
        Xhat(i-1,13),Xhat(i-1,14),Xhat(i-1,15),Xhat(i-1,16),...
        model.I(1,1),model.I(2,2),model.I(3,3),model.m,model.g);
    C = dh;
    % Update step
    [xe, Pe, e, Lambda] = updateStepKF(xn', y(i-1,:)', C, Pn, R, model);
    
    Xhat(i,:) = xe';
    P(:,:,i)  = Pe;
    eta (i,:) = e';
end

% figure
% plot(Xhat(2:end,7:9))
% hold on
% plot(y(1:end,4:6), '--')
%
% figure
% plot(Xhat(2:end,10:12))
% hold on
% plot(y(1:end,7:9), '--')

figure
subplot(311)
shadedErrorBar(t,Xhat(2:end,1),squeeze(2*sqrt(P(1,1,2:end)))','b');
hold on
plot(t, x(:,1), '--b')

subplot(312)
shadedErrorBar(t,Xhat(2:end,2),squeeze(2*sqrt(P(2,2,2:end)))','g');
hold on
plot(t, x(:,2), '--g')

subplot(313)
shadedErrorBar(t,Xhat(2:end,3),squeeze(2*sqrt(P(3,3,2:end)))','r');
hold on
plot(t, x(:,3), '--r')

figure
shadedErrorBar(t,Xhat(2:end,4),squeeze(2*sqrt(P(4,4,2:end)))','b');
hold on
shadedErrorBar(t,Xhat(2:end,5),squeeze(2*sqrt(P(5,5,2:end)))','g');
shadedErrorBar(t,Xhat(2:end,6),squeeze(2*sqrt(P(6,6,2:end)))','r');
plot(t, x(:,4:6), '--')

figure
subplot(221)
shadedErrorBar(t,Xhat(2:end,13),squeeze(2*sqrt(P(13,13,2:end)))','b');
hold on
plot(t, x(:,13), '--b')

subplot(222)
shadedErrorBar(t,Xhat(2:end,14),squeeze(2*sqrt(P(14,14,2:end)))','g');
hold on
plot(t, x(:,14), '--g')

subplot(223)
shadedErrorBar(t,Xhat(2:end,15),squeeze(2*sqrt(P(15,15,2:end)))','r');
hold on
plot(t, x(:,15), '--r')

subplot(224)
shadedErrorBar(t,Xhat(2:end,16),squeeze(2*sqrt(P(16,16,2:end)))','c');
hold on
plot(t, x(:,16), '--c')

for i = 1 : length(t)
    z   (i,:) = q2dcm(   x(i, 13:16))*[0;0;1];
end
tr_param{2} = 10;
for i = 1 : length(t)+1
    % zhat(i,:) = q2dcm(Xhat(i, 13:16))*[0;0;1];
    [zhat(i,:),Pz(:,:,i)] = ut_transform(Xhat(i, 13:16)',P(13:16,13:16,i),@computeVertical, [0;0;1], tr_param);
end

figure
subplot(311)
shadedErrorBar(t,zhat(2:end,1),squeeze(2*sqrt(Pz(1,1,2:end)))','b');
%plot(t,zhat(2:end,1), 'b')
hold on
plot(t,z(:,1), 'b--')
subplot(312)
shadedErrorBar(t,zhat(2:end,2),squeeze(2*sqrt(Pz(2,2,2:end)))','r');
%plot(t,zhat(2:end,2), 'r')
hold on
plot(t,z(:,2), 'r--')
subplot(313)
shadedErrorBar(t,zhat(2:end,3),squeeze(2*sqrt(Pz(3,3,2:end)))','g');
%plot(t,zhat(2:end,3), 'g')
hold on
plot(t,z(:,3), 'g--')
