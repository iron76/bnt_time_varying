% Sestri Leveante, 23 Luglio 2014
%
% This piece of code simulates the problem of estimating dynamic quantities
% for a single rigid body with distributed force/torque measurements and
% distributed gyro/accelerometers measurements. The motion is governed by
% the following differential equation:
%
% m    dv^B    + S(omega^B) (m       v^B) = f^B_1  + ... + f^B_n
%
% I^B domega^B + S(omega^B) (I^B omega^B) = mu^B_1 + ... + mu^B_n
%
% where we defined the following quantities:
%
% I^B    : inertia in the body reference frame
% m      : mass of the rigid body
% f^B_i  : i-th force expressed in the body reference frame
% mu^B_i : i-th torque expressed in the body reference frame
% omega^B: angular velocity expressed in the body reference frame
% v^B    : linear velocity in the body reference frame

clear all
close all
clc

dt      = 0.001;       % sampling time
T       = 1.1;       % time span
sigma_f = 0.01;       % output error
sigma_u = 0.001;      % output error
sigma_a = 0.001;      % output error
n       = 12;         % state dimension
m       = 9;          % output dimension

model.I  = diag([0.05 0.02 0.03]);
model.m  = 20;
model.u  = 1;
model.v  = 0.1;
model.x0 = ones(n,1);
model.dt = dt;

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
% X = [dv^B, domega^b, f^B, mu^B] = [x, f^B, mu^B];
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

Xhat      = zeros(size(x));
Xhat(1,:) = rand(size(model.x0'))*.10 - 5 ;
P(:,:,1)  = eye(n)*100;
Q         = diag([ones(n/2,1).*0.1; ones(n/2,1).*10000]);

Af        = [...
    1/model.m*eye(3), zeros(3,3);
    zeros(3,3)     inv(model.I)];


for i = 2:length(t)+1
    df = rigidBodyDynamicsDerivatives(...
        Xhat(i-1,1) ,Xhat(i-1,2) ,Xhat(i-1,3),...
        Xhat(i-1,4) ,Xhat(i-1,5) ,Xhat(i-1,6),...
        Xhat(i-1,7) ,Xhat(i-1,8) ,Xhat(i-1,9),...
        Xhat(i-1,10),Xhat(i-1,11),Xhat(i-1,12),...        
        model.I(1,1),model.I(2,2),model.I(3,3),model.m);
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
        model.I(1,1),model.I(2,2),model.I(3,3),model.m);
    C = dh;
    % Update step
    [xe, Pe, Lambda] = updateStepKF(xn', y(i-1,:)', C, Pn, R, model);
    
    Xhat(i,:) = xe';
    P(:,:,i)  = Pe;
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
shadedErrorBar(t,Xhat(2:end,1),squeeze(2*sqrt(P(1,1,2:end)))','b');
hold on
shadedErrorBar(t,Xhat(2:end,2),squeeze(2*sqrt(P(2,2,2:end)))','g');
shadedErrorBar(t,Xhat(2:end,3),squeeze(2*sqrt(P(3,3,2:end)))','r');
plot(t, x(:,1:3), '--')

figure
shadedErrorBar(t,Xhat(2:end,4),squeeze(2*sqrt(P(4,4,2:end)))','b');
hold on
shadedErrorBar(t,Xhat(2:end,5),squeeze(2*sqrt(P(5,5,2:end)))','g');
shadedErrorBar(t,Xhat(2:end,6),squeeze(2*sqrt(P(6,6,2:end)))','r');
plot(t, x(:,4:6), '--')



