function [ t,x ] = integrateForwardDynamics( time_span, x0, model )

[n, ~]        = size(x0);
M             = zeros(n,n);
M(1:3 , 1:3)  = eye(3).*model.m;
M(4:6 , 4:6)  = model.I;
M(7:12, 7:12) = eye(6);
odeSettings = odeset('Mass', M);
[t, x]=ode45(@(t, x) rigidBodyDifferentialEquationImplicit(t, x, model), time_span, x0, odeSettings);

end

