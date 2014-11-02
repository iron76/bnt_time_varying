function [h, Sh, index, R] = computeHidden(q, dq, y, o0, p0, do0)

md2q{1} = [];   md2q{1} = [];      md2q{3} = [];   md2q{4} = y(1);
mf{1}   = [];   mf{2}   = y(2:4);  mf{3}   = [];     mf{4}   = [];     mf{5}   = y(8:10);
mu{1}   = [];   mu{2}   = y(5:7);  mu{3}   = [];     mu{4}   = [];     mu{5}   = y(11:13);
mdo{1}  = [];   mdo{2}  = [];      mdo{3}  = [];     mo{4}  = y(14:16)';
mp{1}   = [];   mp{2}   = [];      mp{3}   = [];     mp{4}   = y(17:19)';
fbb     = y(20:22)';
ubb     = y(23:25)';

%dim initialization
d = 3;   %vector dim
n = 4;   %number of links

% number of nodes
k = (n) + (n) +   n   +  n  + ...    
...  |p|    |do|    |d2t|   |c|   
+ (n) + (n) +   n    +   n + ...
...%  |mp|   |mdo|   |md2t|    |mc|
(n+1) + (n+1) +   1   + 1     +   1   + 1 +...
... |u|    |f|     |fbb|   |ubb|   |fee| |uee|
(n+1) + (n+1)     +  n    +   n ;
... |mu|    |mf|    |o|      |mo|

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
Sd2q4 = sd2t;

Sf1 = sF.*eye(d);   Su1 = su.*eye(d);
Sf2 = sF.*eye(d);   Su2 = su.*eye(d);
Sf3 = sF.*eye(d);   Su3 = su.*eye(d);
Sf4 = sF.*eye(d);   Su4 = su.*eye(d);
Sf5 = sF.*eye(d);   Su5 = su.*eye(d);


Sdo1 = sdo.*eye(d); Sp1  = sp.*eye(d); So1  = so.*eye(d);
Sdo2 = sdo.*eye(d); Sp2  = sp.*eye(d); So2  = so.*eye(d);
Sdo3 = sdo.*eye(d); Sp3  = sp.*eye(d); So3  = so.*eye(d);
Sdo4 = sdo.*eye(d); Sp4  = sp.*eye(d); So4  = so.*eye(d);

Sc1  = sp.*eye(d);  Sc2  = sp.*eye(d);
Sc3  = sp.*eye(d);  Sc4  = sp.*eye(d);

sModel   = 1e-6;
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


% indeces
i_p    = [1 2 3 4];        i_do  = [5 6 7 8];
i_d2q  = [9 10 11 12];
i_c    = [13 14 15 16];

m_p    = [17 18 19 20];    m_do  = [21 22 23 24];
m_d2q  = [25 26 27 28];    m_c   = [29 30 31 32];

i_f    = [33 34 35 36 37]; i_u   = [38 39 40 41 42];
i_fbb  = 43;               i_ubb = 44;
i_fee  = 45;               i_uee = 46;
m_f    = [47 48 49 50 51]; m_u   = [52 53 54 55 56];

i_o    = [57 58 59 60];
m_o    = [61 62 63 64];

% sizes
s_p    = [3 3 3 3];        s_do  = [3 3 3 3];
s_d2q  = [1 1 1 1];
s_c    = [3 3 3 3];
s_f    = [3 3 3 3 3];      s_u   = [3 3 3 3 3];
s_fbb  = 3;                s_ubb = 3;
s_fee  = 3;                s_uee = 3;
s_o    = [3 3 3 3];       

% nodes dimesions
ns = [s_p s_do s_d2q s_c s_p s_do s_d2q s_c s_f s_u s_fbb s_ubb s_fee s_uee s_f s_u s_o s_o]; % all 3 dimensional vector for all but the d2q
dnodes = [];                      % no discrete nodes

% conditional relationship
dag = zeros(k);   
dag(i_do(1), i_do(2)) = 1;   
dag(i_do(2), i_do(3)) = 1;   dag(i_d2q(2), i_do(3)) = 1;
dag(i_do(3), i_do(4)) = 1;   dag(i_d2q(3), i_do(4)) = 1;

                              dag(i_do(1), i_p(1))   = 1;   dag(i_d2q(1), i_p(1))  = 1;
dag(i_p(1), i_p(2))    = 1;   dag(i_do(2), i_p(2))   = 1;   dag(i_d2q(2), i_p(2))   = 1;
dag(i_p(2), i_p(3))    = 1;   dag(i_do(3), i_p(3))   = 1;   
dag(i_p(3), i_p(4))    = 1;   dag(i_do(4), i_p(4))   = 1;   

dag(i_p(1), i_c(1))   = 1;   dag(i_do(1), i_c(1))   = 1;
dag(i_p(2), i_c(2))   = 1;   dag(i_do(2), i_c(2))   = 1;
dag(i_p(3), i_c(3))   = 1;   dag(i_do(3), i_c(3))   = 1;
dag(i_p(4), i_c(4))   = 1;   dag(i_do(4), i_c(4))   = 1;


dag(i_p(1), m_p(1) ) = 1;   dag(i_do(1), m_do(1) ) = 1;
dag(i_p(2), m_p(2) ) = 1;   dag(i_do(2), m_do(2) ) = 1;
dag(i_p(3), m_p(3) ) = 1;   dag(i_do(3), m_do(3) ) = 1;
dag(i_p(4), m_p(4) ) = 1;   dag(i_do(4), m_do(4) ) = 1;

dag(i_d2q(1),m_d2q(1)) = 1; dag(i_d2q(2),m_d2q(2)) = 1;
dag(i_d2q(3),m_d2q(3)) = 1; dag(i_d2q(4),m_d2q(4)) = 1;

dag(i_c(1), m_c(1) ) = 1;   dag(i_c(2) , m_c(2)  ) = 1;
dag(i_c(3), m_c(3) ) = 1;   dag(i_c(4) , m_c(4)  ) = 1;


dag(i_fee,  i_f(5))   = 1;   dag(i_uee,  i_u(5))    = 1;

dag(i_f(5), i_f(4))   = 1;   dag(i_c(4), i_f(4))    = 1;
dag(i_f(4), i_f(3))   = 1;   dag(i_c(3), i_f(3))    = 1;   dag(i_fbb, i_f(3))    = 1;
dag(i_f(3), i_f(2))   = 1;   dag(i_c(2), i_f(2))    = 1;
dag(i_f(2), i_f(1))   = 1;   dag(i_c(1), i_f(1))    = 1;

dag(i_do(4), i_u(4))   = 1;   dag(i_f(4), i_u(4))    = 1;
dag(i_f(5),  i_u(4))  = 1;    dag(i_u(5) , i_u(4))   = 1;   

dag(i_do(3), i_u(3))   = 1;   dag(i_f(3), i_u(3))    = 1;  dag(i_ubb, i_u(3))    = 1;
dag(i_f(4),  i_u(3))  = 1;    dag(i_u(4) , i_u(3))   = 1;  dag(i_fbb, i_u(3))    = 1;

dag(i_do(2), i_u(2))   = 1;   dag(i_f(2), i_u(2))    = 1;
dag(i_f(3),  i_u(2))  = 1;    dag(i_u(3) , i_u(2))   = 1;   

dag(i_do(1), i_u(1))  = 1;   dag(i_f(1), i_u(1))    = 1;
dag(i_f(2),  i_u(1))  = 1;   dag(i_u(2),  i_u(1))  = 1;

dag(i_f(1), m_f(1) ) = 1;   dag(i_u(1) , m_u(1)  ) = 1;
dag(i_f(2), m_f(2) ) = 1;   dag(i_u(2) , m_u(2)  ) = 1;
dag(i_f(3), m_f(3) ) = 1;   dag(i_u(3) , m_u(3)  ) = 1;
dag(i_f(4), m_f(4) ) = 1;   dag(i_u(4) , m_u(4)  ) = 1;
dag(i_f(5), m_f(5) ) = 1;   dag(i_u(5) , m_u(5)  ) = 1;

dag(i_o(1), m_o(1) ) = 1;   
dag(i_o(2), m_o(2) ) = 1;   dag(i_o(1), i_o(2)   ) = 1;
dag(i_o(3), m_o(3) ) = 1;   dag(i_o(2), i_o(3)   ) = 1;
dag(i_o(4), m_o(4) ) = 1;   dag(i_o(3), i_o(4)   ) = 1;

% graph init
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

%graph topological reorder such that i-ancestors <  i-descendants

i_p(1)  = find(bnet.order == i_p(1));    i_p(2) = find(bnet.order == i_p(2));
i_p(3)  = find(bnet.order == i_p(3));    i_p(4) = find(bnet.order == i_p(4));
i_do(1) = find(bnet.order == i_do(1));  i_do(2) = find(bnet.order == i_do(2));
i_do(3) = find(bnet.order == i_do(3));  i_do(4) = find(bnet.order == i_do(4));
i_d2q(1)= find(bnet.order == i_d2q(1));i_d2q(2) = find(bnet.order == i_d2q(2));
i_d2q(3)= find(bnet.order == i_d2q(3));i_d2q(4) = find(bnet.order == i_d2q(4));
i_c(1)  = find(bnet.order == i_c(1));    i_c(2) = find(bnet.order == i_c(2));
i_c(3)  = find(bnet.order == i_c(3));    i_c(4) = find(bnet.order == i_c(4));

i_f(1)  = find(bnet.order == i_f(1));    i_f(2) = find(bnet.order == i_f(2));    i_f(3) = find(bnet.order == i_f(3));
i_f(4)  = find(bnet.order == i_f(4));    i_f(5) = find(bnet.order == i_f(5));
i_u(1)  = find(bnet.order == i_u(1));    i_u(2) = find(bnet.order == i_u(2));    i_u(3) = find(bnet.order == i_u(3));
i_u(4)  = find(bnet.order == i_u(4));    i_u(5) = find(bnet.order == i_u(5));    
i_fbb   = find(bnet.order == i_fbb);     i_ubb  = find(bnet.order == i_ubb);    
i_fee   = find(bnet.order == i_fee);     i_uee  = find(bnet.order == i_uee);

m_p(1)  = find(bnet.order == m_p(1));    m_p(2) = find(bnet.order == m_p(2));
m_p(3)  = find(bnet.order == m_p(3));    m_p(4) = find(bnet.order == m_p(4));
m_do(1) = find(bnet.order == m_do(1));  m_do(2) = find(bnet.order == m_do(2));
m_do(3) = find(bnet.order == m_do(3));  m_do(4) = find(bnet.order == m_do(4));
m_d2q(1)= find(bnet.order == m_d2q(1));m_d2q(2) = find(bnet.order == m_d2q(2));
m_d2q(3)= find(bnet.order == m_d2q(3));m_d2q(4) = find(bnet.order == m_d2q(4));
m_c(1)  = find(bnet.order == m_c(1));    m_c(2) = find(bnet.order == m_c(2));
m_c(3)  = find(bnet.order == m_c(3));    m_c(4) = find(bnet.order == m_c(4));

m_f(1)  = find(bnet.order == m_f(1));    m_f(2) = find(bnet.order == m_f(2));    m_f(3) = find(bnet.order == m_f(3));
m_f(4)  = find(bnet.order == m_f(4));    m_f(5) = find(bnet.order == m_f(5));
m_u(1)  = find(bnet.order == m_u(1));    m_u(2) = find(bnet.order == m_u(2));    m_u(3) = find(bnet.order == m_u(3));
m_u(4)  = find(bnet.order == m_u(4));    m_u(5) = find(bnet.order == m_u(5));

i_o(1)  = find(bnet.order == i_o(1));    i_o(2) = find(bnet.order == i_o(2));
i_o(3)  = find(bnet.order == i_o(3));    i_o(4) = find(bnet.order == i_o(4));
m_o(1)  = find(bnet.order == m_o(1));    m_o(2) = find(bnet.order == m_o(2));
m_o(3)  = find(bnet.order == m_o(3));    m_o(4) = find(bnet.order == m_o(4));


identity = eye(k);
new_dag = dag*identity(:,bnet.order);
new_dag = identity(bnet.order, :)*new_dag;
new_ns = ns(bnet.order);

bnet = mk_bnet(new_dag, new_ns, 'discrete', dnodes);
index.p     = i_p;       index.mp  = m_p;
index.do    = i_do;      index.mdo = m_do;
index.d2q   = i_d2q;     index.c   = i_c;
index.md2q  = m_d2q;     index.mc  = m_c;
index.f     = i_f;       index.u   = i_u;
index.fbb   = i_fbb;     index.ubb = i_ubb;
index.fee   = i_fee;     index.uee = i_uee;
index.mf    = m_f;       index.mu  = m_u;
index.o     = i_o;       index.mo  = m_o;

%init distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{index.mdo(1)} = gaussian_CPD(bnet, index.mdo(1),  'mean', zeros(d,1), 'cov', Sdo1, 'clamp_mean', 1, 'clamp_weights', 1,                                    'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mdo(2)} = gaussian_CPD(bnet, index.mdo(2),  'mean', zeros(d,1), 'cov', Sdo2, 'clamp_mean', 1, 'clamp_weights', 1,                                    'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mdo(3)} = gaussian_CPD(bnet, index.mdo(3),  'mean', zeros(d,1), 'cov', Sdo3, 'clamp_mean', 1, 'clamp_weights', 1,                                    'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mdo(4)} = gaussian_CPD(bnet, index.mdo(4),  'mean', zeros(d,1), 'cov', Sdo4, 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight,                'weights', eye(d));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{index.mo(1)}  = gaussian_CPD(bnet, index.mo(1),   'mean', zeros(d,1), 'cov', So1, 'clamp_mean', 1, 'clamp_weights', 1,                                    'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mo(2)}  = gaussian_CPD(bnet, index.mo(2),   'mean', zeros(d,1), 'cov', So2, 'clamp_mean', 1, 'clamp_weights', 1,                                    'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mo(3)}  = gaussian_CPD(bnet, index.mo(3),   'mean', zeros(d,1), 'cov', So3, 'clamp_mean', 1, 'clamp_weights', 1,                                    'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mo(4)}  = gaussian_CPD(bnet, index.mo(4),   'mean', zeros(d,1), 'cov', So4, 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight,                'weights', eye(d));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{index.md2q(1)} = gaussian_CPD(bnet, index.md2q(1),'mean', zeros(1,1), 'cov', Sd2q1, 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight,               'weights', eye(1));
bnet.CPD{index.md2q(2)} = gaussian_CPD(bnet, index.md2q(2),'mean', zeros(1,1), 'cov', Sd2q2, 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight,               'weights', eye(1));
bnet.CPD{index.md2q(3)} = gaussian_CPD(bnet, index.md2q(3),'mean', zeros(1,1), 'cov', Sd2q3, 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight,               'weights', eye(1));
bnet.CPD{index.md2q(4)} = gaussian_CPD(bnet, index.md2q(4),'mean', zeros(1,1), 'cov', Sd2q4, 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight,               'weights', eye(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{index.mf(1)}  = gaussian_CPD(bnet, index.mf(1),   'mean', zeros(d,1), 'cov', Sf1,   'clamp_mean', 1, 'clamp_weights', 1,                                   'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mf(2)}  = gaussian_CPD(bnet, index.mf(2),   'mean', zeros(d,1), 'cov', Sf2,   'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight,               'weights', eye(d));
bnet.CPD{index.mf(3)}  = gaussian_CPD(bnet, index.mf(3),   'mean', zeros(d,1), 'cov', Sf3,   'clamp_mean', 1, 'clamp_weights', 1,                                   'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mf(4)}  = gaussian_CPD(bnet, index.mf(4),   'mean', zeros(d,1), 'cov', Sf4,   'clamp_mean', 1, 'clamp_weights', 1,                                   'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mf(5)}  = gaussian_CPD(bnet, index.mf(5),   'mean', zeros(d,1), 'cov', Sf5,   'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight,               'weights', eye(d));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{index.mu(1)}  = gaussian_CPD(bnet, index.mu(1),   'mean', zeros(d,1), 'cov', Su1, 'clamp_mean', 1, 'clamp_weights', 1,                                     'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mu(2)}  = gaussian_CPD(bnet, index.mu(2),   'mean', zeros(d,1), 'cov', Su2, 'clamp_mean', 1, 'clamp_weights', 1,   'cov_prior_weight', covPriorWeight,               'weights', eye(d));
bnet.CPD{index.mu(3)}  = gaussian_CPD(bnet, index.mu(3),   'mean', zeros(d,1), 'cov', Su3, 'clamp_mean', 1, 'clamp_weights', 1,                                     'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mu(4)}  = gaussian_CPD(bnet, index.mu(4),   'mean', zeros(d,1), 'cov', Su4, 'clamp_mean', 1, 'clamp_weights', 1,                                     'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mu(5)}  = gaussian_CPD(bnet, index.mu(5),   'mean', zeros(d,1), 'cov', Su5, 'clamp_mean', 1, 'clamp_weights', 1,   'cov_prior_weight', covPriorWeight,               'weights', eye(d));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{index.mp(1)}  = gaussian_CPD(bnet, index.mp(1),   'mean', zeros(d,1), 'cov', Sp1, 'clamp_mean', 1, 'clamp_weights', 1,                                     'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mp(2)}  = gaussian_CPD(bnet, index.mp(2),   'mean', zeros(d,1), 'cov', Sp2, 'clamp_mean', 1, 'clamp_weights', 1,                                     'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mp(3)}  = gaussian_CPD(bnet, index.mp(3),   'mean', zeros(d,1), 'cov', Sp3, 'clamp_mean', 1, 'clamp_weights', 1,                                     'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mp(4)}  = gaussian_CPD(bnet, index.mp(4),   'mean', zeros(d,1), 'cov', Sp4, 'clamp_mean', 1, 'clamp_weights', 1,   'cov_prior_weight', covPriorWeight,               'weights', eye(d));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{index.mc(1)}  = gaussian_CPD(bnet, index.mc(1),   'mean', zeros(d,1), 'cov', Sc1, 'clamp_mean', 1, 'clamp_weights', 1,                                     'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mc(2)}  = gaussian_CPD(bnet, index.mc(2),   'mean', zeros(d,1), 'cov', Sc2, 'clamp_mean', 1, 'clamp_weights', 1,                                     'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mc(3)}  = gaussian_CPD(bnet, index.mc(3),   'mean', zeros(d,1), 'cov', Sc3, 'clamp_mean', 1, 'clamp_weights', 1,                                     'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mc(4)}  = gaussian_CPD(bnet, index.mc(4),   'mean', zeros(d,1), 'cov', Sc4, 'clamp_mean', 1, 'clamp_weights', 1,                                     'clamp_cov', 1, 'weights', eye(d));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{index.fbb}    = gaussian_CPD(bnet, index.fbb,     'mean', zeros(3,1), 'cov', Sfbb,'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
bnet.CPD{index.ubb}    = gaussian_CPD(bnet, index.ubb,     'mean', zeros(3,1), 'cov', Subb,'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{index.fee}    = gaussian_CPD(bnet, index.fee,     'mean', zeros(3,1), 'cov', Sfee,'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
bnet.CPD{index.uee}    = gaussian_CPD(bnet, index.uee,     'mean', zeros(3,1), 'cov', Suee,'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);


%unknowns
bnet.CPD{index.d2q(1)}   = gaussian_CPD(bnet, index.d2q(1),'mean', zeros(1,1), 'cov', sUnknown.*eye(1), 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
bnet.CPD{index.d2q(2)}   = gaussian_CPD(bnet, index.d2q(2),'mean', zeros(1,1), 'cov', sUnknown.*eye(1), 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
bnet.CPD{index.d2q(3)}   = gaussian_CPD(bnet, index.d2q(3),'mean', zeros(1,1), 'cov', sUnknown.*eye(1), 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
bnet.CPD{index.d2q(4)}   = gaussian_CPD(bnet, index.d2q(4),'mean', zeros(1,1), 'cov', sUnknown.*eye(1), 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);

ns = new_ns;

load columnImplicit.mat

c1 = cos(q(3));  c2 = cos(q(4));
s1 = sin(q(3));  s2 = sin(q(4));

% R23
R(1:3,1:3,3) = [s1 c1 0;
      -c1  s1 0;
      0    0 1];
% R34
R(1:3,1:3,4) = [c2 -s2 0;
      s2  c2 0;
      0    0 1];

R(1:3,1:3,5) = [1 0 0;
      0 1 0;
      0 0 1];  
% R01
R(1:3,1:3,1) = [1  0 0;
                0  0 1;
                0 -1 0];
% R12  
R(1:3,1:3,2) = [0 0 -1;
                0 1  0;
                1 0  0];

dt  = [dq(1), dq(2),  dq(3),  dq(4)];            

%Useful quantities for implicit NE
param = [l1 l2 m1 m2 r1 r2 g I1z I2z];
xK = [ q'; dq'];
Ci  = CY(xK, param); 
bi  = bY(xK, param);

indC = 1;
%Kinematic recursion
for i = 1:4
    
    j = (1+d*(i-1)):(d*i);
    j_prev = (1+d*(i-2)):(d*(i-1));
    
    Ri = R(:,:,i);
    if i ~= 1
        index_do_prev = index.do(i-1);
        index_p_prev  =  index.p(i-1);
    end
    
    if i == 3 || i == 4
        
        
        %% o(:,i)  = Ri'*(o_prev+dq(i)*z0);
        W = Ci(indC:indC+2, colum.o(j_prev));
        b = bi(indC:indC+2);
        bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;
        
        %% do(:,i) = Ri'*(do_prev+d2q(i)*z0+dq(i)*cross(o_prev,z0));
        b  = bi(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.d2t(i));
        W2 = Ci(indC:indC+2, colum.do(j_prev));
        W = commuteOnOrder2(W1, W2, index.d2q(i), index_do_prev);
        bnet.CPD{index.do(i)}  = gaussian_CPD(bnet, index.do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;

        %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
        b  = bi(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.p(j_prev));
        W2 = Ci(indC:indC+2, colum.do(j));
        W = commuteOnOrder2(W1, W2, index_p_prev, index.do(i));
        bnet.CPD{index.p(i)}  = gaussian_CPD(bnet, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;
        
        
        %% c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));
        b  = bi(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.p(j)) ;
        W2 = Ci(indC:indC+2, colum.do(j));
        
        W = commuteOnOrder2(W1, W2, index.p(i), index.do(i));
        bnet.CPD{index.c(i)}  = gaussian_CPD(bnet, index.c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;

    else
        
        %% o(:,i)  = Ri'*o_prev;
        if i ~= 1
            W = Ci(indC:indC+2, colum.o(j_prev));
            b = bi(indC:indC+2);
            bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        else
            b = bi(indC:indC+2);
            bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', b, 'cov', sModel.*eye(d));
        end
        indC = indC+3;
        
        %% do(:,i) = Ri'*do_prev;
        if i == 1
            b = bi(indC:indC+2);
            bnet.CPD{index.do(i)}  = gaussian_CPD(bnet, index.do(i),  'mean', b, 'cov', sModel.*eye(d));
        else
            W = Ci(indC:indC+2, colum.do(j_prev));
            b = bi(indC:indC+2);
            bnet.CPD{index.do(i)}  = gaussian_CPD(bnet, index.do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        end
        indC = indC+3;
        
        %% p(:,i)  = Ri'*(p_prev+d2t(i)*z0)+2*dt(i)*cross(o(:,i), Ri'*z0)+cross(do(:,i),ri0_i(:,i))+cross(o(:,i),cross(o(:,i),ri0_i(:,i)));
        b  = bi(indC:indC+2);
        if i == 1
            W1 = Ci(indC:indC+2, colum.d2t(i));
            W2 = Ci(indC:indC+2, colum.do(j));
            W = commuteOnOrder2(W1, W2, index.d2q(i), index.do(i));
        else
            W1 = Ci(indC:indC+2, colum.p(j_prev));
            W2 = Ci(indC:indC+2, colum.d2t(i));
            W3 = Ci(indC:indC+2, colum.do(j));
            W = commuteOnOrder3(W1, W2, W3, index_p_prev, index.d2q(i), index.do(i));
        end
        bnet.CPD{index.p(i)}  = gaussian_CPD(bnet, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;
        
        
                
        %% c(:,i)  = p(:,i) + cross(do(:,i), riC_i(:,i)) +cross(o(:,i),cross(o(:,i),riC_i(:,i)));
        b = bi(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.p(j));
        W2 = Ci(indC:indC+2, colum.do(j));
        W = commuteOnOrder2(W1, W2, index.p(i), index.do(i));
        bnet.CPD{index.c(i)}  = gaussian_CPD(bnet, index.c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W); 
        indC = indC+3;

    end
end

j = (1+d*4):(d*5);

% f(5) = Ree * fee;
b  = bi(indC:indC+2);
Ree = Ci(indC:indC+2, colum.fee);
bnet.CPD{index.f(5)}  = gaussian_CPD(bnet, index.f(5),  'mean', b, 'cov', sModel.*eye(d), 'weights', Ree);
indC = indC+3;

% u(5) = Ree * uee;
b  = bi(indC:indC+2);
Ree = Ci(indC:indC+2, colum.uee);
bnet.CPD{index.u(5)}  = gaussian_CPD(bnet, index.u(5),  'mean', b, 'cov', sModel.*eye(d), 'weights', Ree);
indC = indC+3;


%Dynamic recursion
for i = 4:-1:1

    j = (1+d*(i-1)):(d*i);
    j_post = (1+d*(i)):(d*(i+1));

    Ri = R(:,:,i+1);
    Ii = I(:,:,i);
    if i == 3
        %% f(:,i)  = Ri*f(:,i+1) + m(i,1).*c(:,i) + Rbb*fbb;
        %  fbbNE   = Rbb * fbb;
        %  ubbNE   = Rbb * ubb;
        %  Rbb     = -1.*(Rb0*R(:,:,1)*R(:,:,2)*R(:,:,3))';        
        b  = bi(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.f(j_post));
        W2 = Ci(indC:indC+2, colum.c(j));
        W3 = Ci(indC:indC+2, colum.fbb);
        
        W = commuteOnOrder3(W1, W2, W3, index.f(i+1), index.c(i), index.fbb);
        bnet.CPD{index.f(i)}  = gaussian_CPD(bnet, index.f(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;
        
        
        %% u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
        %    cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i)) + cross(-rc(:,i)-r0(:,i), Rbb*fbb) + Rbb*ubb;
        b = bi(indC:indC+2);
        Wt(:, :, 1) = Ci(indC:indC+2, colum.do(j));
        Wt(:, :, 2) = Ci(indC:indC+2, colum.f(j));
        Wt(:, :, 3) = Ci(indC:indC+2, colum.f(j_post));
        Wt(:, :, 4) = Ci(indC:indC+2, colum.u(j_post));
        Wt(:, :, 5) = Ci(indC:indC+2, colum.fbb);
        Wt(:, :, 6) = Ci(indC:indC+2, colum.ubb);
        [~,ind] = sort([index.do(i), index.f(i), index.f(i+1),  index.u(i+1), index.fbb, index.ubb]);
        W = [Wt(:,:,ind(1)) Wt(:,:,ind(2)) Wt(:,:,ind(3)) Wt(:,:,ind(4)) Wt(:,:,ind(5)) Wt(:,:,ind(6))];
        
        bnet.CPD{index.u(i)}  = gaussian_CPD(bnet, index.u(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);                
        indC = indC+3;

        
    else
        %% f(:,i)  = Ri*f(:,i+1) + m(i,1).*c(:,i);
        b =  bi(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.f(j_post));
        W2 = Ci(indC:indC+2, colum.c(j));
        W = commuteOnOrder2(W1, W2, index.f(i+1), index.c(i));
        bnet.CPD{index.f(i)}  = gaussian_CPD(bnet, index.f(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;
        
        %% u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
        %    cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i))
        b = bi(indC:indC+2);
        Wt(:, :, 1) = Ci(indC:indC+2, colum.do(j));
        Wt(:, :, 2) = Ci(indC:indC+2, colum.f(j));
        Wt(:, :, 3) = Ci(indC:indC+2, colum.f(j_post));
        Wt(:, :, 4) = Ci(indC:indC+2, colum.u(j_post));
        [~,ind] = sort([index.do(i), index.f(i), index.f(i+1),  index.u(i+1)]);
        W = [Wt(:,:,ind(1)) Wt(:,:,ind(2)) Wt(:,:,ind(3)) Wt(:,:,ind(4))];
        
        bnet.CPD{index.u(i)}  = gaussian_CPD(bnet, index.u(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);                
        indC = indC+3;
    end
end


engine = jtree_inf_engine(bnet);

% set conditional nodes
evidence = cell(1,k);

% evidence{index.md2q(3)}  = md2q{3};
evidence{index.md2q(4)}  = md2q{4};
evidence{index.mf(2)}    = R(1:3,1:3,3)*mf{2};
evidence{index.mu(2)}    = R(1:3,1:3,3)*mu{2};
evidence{index.mf(5)}    = mf{5};
evidence{index.mu(5)}    = mu{5};
evidence{index.mo(4)}    = mo{4};
evidence{index.mp(4)}    = mp{4};
evidence{index.fbb}      = fbb;
evidence{index.ubb}      = ubb;

[engine, ll] = enter_evidence(engine, evidence);

for i = 1 : n
    marg = marginal_nodes(engine, index.p(i));
    h{index.p(i)}  = marg.mu;
    Sh{index.p(i)} = marg.Sigma;
end

for i = 1 : n
    marg = marginal_nodes(engine, index.o(i));
    h{index.o(i)}  = marg.mu;
    Sh{index.o(i)} = marg.Sigma;
end

for i = 1 : n
    marg = marginal_nodes(engine, index.do(i));
    h{index.do(i)}  = marg.mu;
    Sh{index.do(i)} = marg.Sigma;
end

for i = 1 : n
    marg = marginal_nodes(engine, index.d2q(i));
    h{index.d2q(i)}  = marg.mu;
    Sh{index.d2q(i)} = marg.Sigma;
end

for i = 1 : n
    marg = marginal_nodes(engine, index.c(i));
    h{index.c(i)}  = marg.mu;
    Sh{index.c(i)} = marg.Sigma;
end

for i = 1 : n+1
    marg = marginal_nodes(engine, index.f(i));
    h{index.f(i)}  = marg.mu;
    Sh{index.f(i)} = marg.Sigma;
end

for i = 1 : n+1
    marg = marginal_nodes(engine, index.u(i));
    h{index.u(i)}  = marg.mu;
    Sh{index.u(i)} = marg.Sigma;
end

marg = marginal_nodes(engine, index.fbb);
h{index.fbb}  = marg.mu;
Sh{index.fbb} = marg.Sigma;

marg = marginal_nodes(engine, index.ubb);
h{index.ubb}  = marg.mu;
Sh{index.ubb} = marg.Sigma;

marg = marginal_nodes(engine, index.fee);
h{index.fee}  = marg.mu;
Sh{index.fee} = marg.Sigma;

marg = marginal_nodes(engine, index.uee);
h{index.uee}  = marg.mu;
Sh{index.uee} = marg.Sigma;
