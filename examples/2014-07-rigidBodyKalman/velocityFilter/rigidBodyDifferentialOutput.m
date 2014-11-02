function [dx, f, mu] = rigidBodyDifferentialOutput(t,x,p)

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

v_B     = x(1:3, 1);
omega_B = x(4:6, 1); 

I_B  = p.I;
m    = p.m;
f_B  = p.f;
mu_B = p.mu;
ts   = p.ts;

f_B_t  = interp1(ts,   f_B , t)';
mu_B_t = interp1(ts,   mu_B, t)';

dv_B     =          -S(omega_B) * v_B + 1/m * f_B_t;
domega_B =   I_B \ (-S(omega_B) * (I_B * omega_B) + mu_B_t);

dx  = [dv_B; domega_B];
