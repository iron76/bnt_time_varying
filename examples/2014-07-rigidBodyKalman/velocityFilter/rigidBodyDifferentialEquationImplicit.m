function f = rigidBodyDifferentialEquationImplicit(t,x,p)

% This piece of code computes the rigid body differential equations. The
% equation is defined as follows:
%
% m    dv^B    = -S(omega^B) (m       v^B) + f^B_1  + ... + f^B_n
% 
% I^B domega^B = -S(omega^B) (I^B omega^B) + mu^B_1 + ... + mu^B_n
% 

v_B     = x(1:3  , 1);
omega_B = x(4:6  , 1); 
f_B_t   = x(7:9  , 1);
mu_B_t  = x(10:12, 1); 


I_B  = p.I;
m    = p.m;
u    = p.u;
v    = p.v;

u_t  = [cos(t); sin(t); sin(2*t)]*u;
v_t  = [sin(t); cos(t); sin(2*t)]*v;

fLin     = - S(omega_B) * (m   * v_B    ) + f_B_t;
fAng     = - S(omega_B) * (I_B * omega_B) + mu_B_t;
df_B     =                                     u_t;
dmu_B    =                                     v_t;

f  = [fLin; fAng; df_B; dmu_B];
