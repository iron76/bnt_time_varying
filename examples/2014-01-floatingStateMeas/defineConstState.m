%dim initialization
d  = 3;   %vector dim
n  = 4;   %number of links
dx = 2*n; %state dimension


% number of nodes
k =  n +  n    +   n   +  n  + ...    
...  |p|    |do|    |d2q|   |c|   
+      n +     n    +   n    +   n + ...
...%  |mp|   |mdo|   |md2q|    |mc|
(n+1) + (n+1) +   1   + 1     +   1   + 1 +...
... |u|    |f|     |fbb|   |ubb|   |fee| |uee|
(n+1) + (n+1)   + n +  n +  1  + 1;
... |mu|    |mf| |o|  |mo| |x|   |mx|

% Rb0 from NE base to x-y plane
Rb0 = [0  0 1;
       0  1 0;
       -1 0 0];

%prior on covariance
covPriorWeight = 0.2;

%uncertainties
sdo      = 0.7805; 
sp       = 0.01;
sd2t     = 10;
sF       = 10;
su       = 0.0031;
so       = 0.01;

Sd2q1 = sd2t;
Sd2q2 = sd2t;
Sd2q3 = sd2t;
Sd2q4 = sd2t;    Sdq4 = 10*0.1*pi/180;       Sq4 = 0.1*pi/180;

Sf1 = sF.*eye(d);   Su1 = su.*eye(d);
Sf2 = sF.*eye(d);   Su2 = su.*eye(d);
Sf3 = sF.*eye(d);   Su3 = su.*eye(d);
Sf4 = sF.*eye(d);   Su4 = su.*eye(d);
Sf5 = sF.*eye(d);   Su5 = su.*eye(d);

Smx = diag([Sq4, Sdq4]);
Wx  = [0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 1]; 

Sdo1 = sdo.*eye(d); Sp1  = sp.*eye(d); So1  = so.*eye(d);
Sdo2 = sdo.*eye(d); Sp2  = sp.*eye(d); So2  = so.*eye(d);
Sdo3 = sdo.*eye(d); Sp3  = sp.*eye(d); So3  = so.*eye(d);
Sdo4 = sdo.*eye(d); Sp4  = sp.*eye(d); So4  = so.*eye(d);

Sc1  = sp.*eye(d);  Sc2  = sp.*eye(d);
Sc3  = sp.*eye(d);  Sc4  = sp.*eye(d);

sModel   = 5e-3;
sUnknown = 1e6;

Sdo0 = sModel.*eye(d); Sp0  = sModel.*eye(d);
Sfbb = sModel.*eye(d); Subb = sModel.*eye(d);
Sfee = sModel.*eye(d); Suee = sModel.*eye(d);

load cov_hat.mat

%constant parameters
l1 = 0.2236;     l2 = 0.213;
lc1 = -0.5*l1;   lc2 = -0.5*l2;
g = 9.81;
m1 = 0.754+0.526+2.175;     m2 = 1.264+0.746+0.010;
I1z = 0.00001;   I2z = 0.00005;

g0 = [0; g; 0];
z0 = [0; 0; 1];
%Dynamic parameters
I(1:3,1:3,1) = [0 0 0; 0 0 0; 0 0 0];
I(1:3,1:3,2) = [0 0 0; 0 0 0; 0 0 0];
I(1:3,1:3,3) = [0 0 0; 0 0 0; 0 0 I1z];
I(1:3,1:3,4) = [0 0 0; 0 0 0; 0 0 I2z];
m  = [0; 0; m1; m2];

%Kinmeatic parameters
r1c = [lc1; 0; 0];
r2c = [lc2; 0; 0];
r10 = [l1 ; 0; 0];
r20 = [l2 ; 0; 0];
r0  = [zeros(3,1), zeros(3,1), r10,  r20];
rc  = [zeros(3,1), zeros(3,1), r1c,  r2c];
r1  = lc1 + l1;
r2  = lc2 + l2;

%FT sensor rotation
R = Rz(pi/2)*Rx(pi);
Adj_IMU = [R zeros(3,3); zeros(3,3) R];

%IMU sensor rotation
R = Rz(pi)*Rx(pi/2);
Adj_FT = [R zeros(3,3); zeros(3,3) R];


%Exec options
REMOVE_OFFSETS = 1;

%initial position offset in degrees
q0_offset = [90 0 0 0 0 0]';    
