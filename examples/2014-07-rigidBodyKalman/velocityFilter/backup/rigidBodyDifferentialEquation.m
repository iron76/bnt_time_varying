function dx = rigidBodyDifferentialEquation(t,x,p)

% This piece of code computes the rigid body differential equations. The
% equation is defined as follows:
%
% m    dv^B    + S(omega^B) (m       v^B) = f^B_1  + ... + f^B_n
% 
% I^B domega^B + S(omega^B) (I^B omega^B) = mu^B_1 + ... + mu^B_n
% 
% and therefore:
%
%      dv^B = - S(omega^B) (        v^B) + 1/m f^B_1  + ... + 1/m f^B_n
% 
%  domega^B = - inv(I^B) (S(omega^B) (I^B omega^B) + mu^B_1 + ... + mu^B_n)
%
% Since also forces and torques will be part of the system state for the
% Kalman filter, we add also a differential equation for the time evolution
% of force and torques. 
%
%    df^B  = u
%
%    dmu^B = v


v_B     = x(1:3  , 1);
omega_B = x(4:6  , 1);
f_B     = x(7:9  , 1);
mu_B    = x(10:12, 1);


I_B  = p.I;
m    = p.m;
u    = p.u;
v    = p.v;
ts   = p.ts;

u_t  = interp1(ts,   u, t)';
v_t  = interp1(ts,   v, t)';

dv_B     =          -S(omega_B) * v_B + 1/m * f_B;
domega_B =   I_B \ (-S(omega_B) * (I_B * omega_B) + mu_B);
df_B     =                                     u_t;
dmu_B    =                                     v_t;

dx  = [dv_B; domega_B; df_B; dmu_B];
