function M = massMatrix( ~, ~, md )
I_B  = md.I;
m    = md.m;

M    = ...
    [eye(3).*m, zeros(3,3)
    zeros(3,3) I_B];
end

