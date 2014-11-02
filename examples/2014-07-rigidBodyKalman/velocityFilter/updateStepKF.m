function [xe, Pe, Lambda] = updateStepKF(x, y, C, P, R, model)

[n, ~] = size(x);

Lambda = C*P*C'+R;
Mn = P*C'/Lambda;                              % Innovation
xe = x + Mn*(y-rigidBodyOutput([],x,model));   % x[n|n]
Pe = (eye(n)-Mn*C)*P;                          % P[n|n]

end