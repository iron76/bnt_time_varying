%dim initialization
d = 3;   %vector dim
n = 2;   %number of links

% number of nodes
k = n*6+4;

%uncertainties
sdo = 0.7805; 
sp  = 0.01;
sd2t= 10;
sf  = 0.1;
su  = 0.0031;
sModel = .0001;

Sd2t1 = sd2t;
Sd2t2 = sd2t;
Sf1 = sf.*eye(d,d);   Su1 = su.*eye(d,d);
Sf2 = sf.*eye(d,d);   Su2 = su.*eye(d,d);
Sf3 = sf.*eye(d,d);   Su3 = su.*eye(d,d);

Sdo0 = sdo.*eye(d,d); Sp0  = sp.*eye(d,d);
Sdo1 = sdo.*eye(d,d); Sp1  = sp.*eye(d,d);
Sdo2 = sdo.*eye(d,d); Sp2  = sp.*eye(d,d);

%constant parameters
l1 = 0.2236;     l2 = 0.213;
lc1 = -0.5*l1;   lc2 = -0.5*l2;
g = 9.81;
m1 = 0.754+0.526+2.175;     m2 = 1.264+0.746+0.010;
I1z = 0.00001;   I2z = 0.00005;

g0 = [0; g; 0];
z0 = [0; 0; 1];
%Dynamic parameters
I1 = [0 0 0; 0 0 0; 0 0 I1z];
I2 = [0 0 0; 0 0 0; 0 0 I2z];
m = [m1; m2];

%Kinmeatic parameters
r1c = [lc1; 0; 0];
r2c = [lc2; 0; 0];
r10 = [l1 ; 0; 0];
r20 = [l2 ; 0; 0];
r0  = [r10,  r20];
rc  = [r1c,  r2c];

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




