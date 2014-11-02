function [bnet, index] = init2DOFnet3(n, k, d)

% graph initialization. graph contains the following nodes:
%
% - linear and angular link accelerations p, do including base [2*(n+1)]
%     p0, do0, p = [p1, p2]; do = [do1, do2],
% - forces and torques including end effector f, u [2*(n+1)]
%     f = [f1, f2, f3];  u = [u1, u2, u3];
% - joint acceletations d2t1, d2t2 [n]
%     d2t1, d2t2
% - center of mass acceleration [n]
%     c = [c1, c2]

defineConst2

% indeces
i_p0   = 1;       i_do0 = 2;
i_p    = [3 4];   i_do  = [5 6];
i_f    = [7 8 9]; i_u   = [10 11 12];
i_d2t  = [13 14]; i_c   = [15 16];
%error measurements
m_p0   = 17;       m_do0 = 18;
m_p    = [19 20];   m_do  = [21 22];
m_f    = [23 24 25]; m_u   = [26 27 28];
m_d2t  = [29 30]; m_c   = [31 32];

% nodes dimesions
ns = [d*ones(1, 4*n+4), ones(1,n), d*ones(1, n), d*ones(1, 4*n+4), ones(1,n), d*ones(1, n) ]; % all 3 dimensional vector for all but the d2t
dnodes = [];                      % no discrete nodes

% conditional relationship
dag = zeros(k);
dag(i_do0  , i_do(1)) = 1;   dag(i_d2t(1), i_do(1)) = 1;
dag(i_do(1), i_do(2)) = 1;   dag(i_d2t(2), i_do(2)) = 1;

dag(i_p0  , i_p(1))    = 1;   dag(i_do(1), i_p(1))   = 1;
dag(i_p(1), i_p(2))    = 1;   dag(i_do(2), i_p(2))   = 1;

dag(i_p(1), i_c(1))   = 1;   dag(i_do(1), i_c(1))   = 1;
dag(i_p(2), i_c(2))   = 1;   dag(i_do(2), i_c(2))   = 1;

dag(i_f(3), i_f(2))   = 1;   dag(i_c(2), i_f(2))    = 1;
dag(i_f(2), i_f(1))   = 1;   dag(i_c(1), i_f(1))    = 1;

dag(i_do(2), i_u(2))   = 1;   dag(i_f(2), i_u(2))    = 1;
dag(i_f(3),  i_u(2))  = 1;    dag(i_u(3) , i_u(2))   = 1;   

dag(i_do(1), i_u(1))  = 1;   dag(i_f(1), i_u(1))    = 1;
dag(i_f(2),  i_u(1))  = 1;   dag(i_u(2),  i_u(1))  = 1;

dag(i_p0  , m_p0   ) = 1;   dag(i_do0  , m_do0   ) = 1;

dag(i_p(1), m_p(1) ) = 1;   dag(i_do(1), m_do(1) ) = 1;
dag(i_p(2), m_p(2) ) = 1;   dag(i_do(2), m_do(2) ) = 1;

dag(i_f(1), m_f(1) ) = 1;   dag(i_u(1) , m_u(1)  ) = 1;
dag(i_f(2), m_f(2) ) = 1;   dag(i_u(2) , m_u(2)  ) = 1;
dag(i_f(3), m_f(3) ) = 1;   dag(i_u(3) , m_u(3)  ) = 1;

dag(i_d2t(1),m_d2t(1)) = 1; dag(i_d2t(2),m_d2t(2)) = 1;

dag(i_c(1), m_c(1) ) = 1;   dag(i_c(2) , m_c(2)  ) = 1;


% graph init
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

%graph topological reorder such that i-ancestors <  i-descendants

i_p0    = find(bnet.order == i_p0);       i_do0 = find(bnet.order == i_do0);
i_p(1)  = find(bnet.order == i_p(1));    i_p(2) = find(bnet.order == i_p(2));
i_do(1) = find(bnet.order == i_do(1));  i_do(2) = find(bnet.order == i_do(2));
i_f(1)  = find(bnet.order == i_f(1));    i_f(2) = find(bnet.order == i_f(2));    i_f(3) = find(bnet.order == i_f(3));
i_u(1)  = find(bnet.order == i_u(1));    i_u(2) = find(bnet.order == i_u(2));    i_u(3) = find(bnet.order == i_u(3));
i_d2t(1)= find(bnet.order == i_d2t(1));i_d2t(2) = find(bnet.order == i_d2t(2));
i_c(1)  = find(bnet.order == i_c(1));    i_c(2) = find(bnet.order == i_c(2));

m_p0    = find(bnet.order == m_p0);       m_do0 = find(bnet.order == m_do0);
m_p(1)  = find(bnet.order == m_p(1));    m_p(2) = find(bnet.order == m_p(2));
m_do(1) = find(bnet.order == m_do(1));  m_do(2) = find(bnet.order == m_do(2));
m_f(1)  = find(bnet.order == m_f(1));    m_f(2) = find(bnet.order == m_f(2));    m_f(3) = find(bnet.order == m_f(3));
m_u(1)  = find(bnet.order == m_u(1));    m_u(2) = find(bnet.order == m_u(2));    m_u(3) = find(bnet.order == m_u(3));
m_d2t(1)= find(bnet.order == m_d2t(1));m_d2t(2) = find(bnet.order == m_d2t(2));
m_c(1)  = find(bnet.order == m_c(1));    m_c(2) = find(bnet.order == m_c(2));

identity = eye(k);
new_dag = dag*identity(:,bnet.order);
new_dag = identity(bnet.order, :)*new_dag;
new_ns = ns(bnet.order);

bnet = mk_bnet(new_dag, new_ns, 'discrete', dnodes);
index.p0    = i_p0;      index.p  = i_p;
index.do0   = i_do0;     index.do = i_do;
index.f     = i_f;       index.u  = i_u;
index.d2t   = i_d2t;     index.c  = i_c;

index.mp0    = m_p0;      index.mp  = m_p;
index.mdo0   = m_do0;     index.mdo = m_do;
index.mf     = m_f;       index.mu  = m_u;
index.md2t   = m_d2t;     index.mc  = m_c;

%init distributions
bnet.CPD{index.mp0}    = gaussian_CPD(bnet, index.mp0,     'mean', zeros(d,1), 'cov', Sp0,   'weights', eye(d));
bnet.CPD{index.mp(1)}  = gaussian_CPD(bnet, index.mp(1),   'mean', zeros(d,1), 'cov', Sp1,   'weights', eye(d));
bnet.CPD{index.mp(2)}  = gaussian_CPD(bnet, index.mp(2),   'mean', zeros(d,1), 'cov', Sp2,   'weights', eye(d));

bnet.CPD{index.mdo0}   = gaussian_CPD(bnet, index.mdo0,    'mean', zeros(d,1), 'cov', Sdo0,  'weights', eye(d));
bnet.CPD{index.mdo(1)} = gaussian_CPD(bnet, index.mdo(1),  'mean', zeros(d,1), 'cov', Sdo1,  'weights', eye(d));
bnet.CPD{index.mdo(2)} = gaussian_CPD(bnet, index.mdo(2),  'mean', zeros(d,1), 'cov', Sdo2,  'weights', eye(d));

bnet.CPD{index.mf(1)}  = gaussian_CPD(bnet, index.mf(1),   'mean', zeros(d,1), 'cov', Sf1,   'weights', eye(d));
bnet.CPD{index.mf(2)}  = gaussian_CPD(bnet, index.mf(2),   'mean', zeros(d,1), 'cov', Sf2,   'weights', eye(d));
bnet.CPD{index.mf(3)}  = gaussian_CPD(bnet, index.mf(3),   'mean', zeros(d,1), 'cov', Sf3,   'weights', eye(d));

bnet.CPD{index.mu(1)}  = gaussian_CPD(bnet, index.mu(1),   'mean', zeros(d,1), 'cov', Su1,   'weights', eye(d));
bnet.CPD{index.mu(2)}  = gaussian_CPD(bnet, index.mu(2),   'mean', zeros(d,1), 'cov', Su2,   'weights', eye(d));
bnet.CPD{index.mu(3)}  = gaussian_CPD(bnet, index.mu(3),   'mean', zeros(d,1), 'cov', Su3,   'weights', eye(d));

bnet.CPD{index.md2t(1)} = gaussian_CPD(bnet, index.md2t(1),'mean', zeros(1,1), 'cov', Sd2t1, 'weights', eye(1));
bnet.CPD{index.md2t(2)} = gaussian_CPD(bnet, index.md2t(2),'mean', zeros(1,1), 'cov', Sd2t2, 'weights', eye(1));

bnet.CPD{index.mc(1)}  = gaussian_CPD(bnet, index.mc(1),   'mean', zeros(d,1), 'cov', Sc1,   'weights', eye(d));
bnet.CPD{index.mc(2)}  = gaussian_CPD(bnet, index.mc(2),   'mean', zeros(d,1), 'cov', Sc2,   'weights', eye(d));

%unknowns
bnet.CPD{index.do0}      = gaussian_CPD(bnet, index.do0,      'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.p0}       = gaussian_CPD(bnet, index.p0,       'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.d2t(1)}   = gaussian_CPD(bnet, index.d2t(1),   'mean', zeros(1,1),    'cov', sUnknown.*eye(1,1));
bnet.CPD{index.d2t(2)}   = gaussian_CPD(bnet, index.d2t(2),   'mean', zeros(1,1),    'cov', sUnknown.*eye(1,1));
bnet.CPD{index.f(3)}     = gaussian_CPD(bnet, index.f(3),     'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.u(3)}     = gaussian_CPD(bnet, index.u(3),     'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.do(1)}    = gaussian_CPD(bnet, index.do(1),    'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.c(1)}     = gaussian_CPD(bnet, index.c(1),     'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.c(2)}     = gaussian_CPD(bnet, index.c(2),     'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.p(1)}     = gaussian_CPD(bnet, index.p(1),     'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.p(2)}     = gaussian_CPD(bnet, index.p(2),     'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.f(1)}     = gaussian_CPD(bnet, index.f(1),     'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.u(1)}     = gaussian_CPD(bnet, index.u(1),     'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.do(2)}    = gaussian_CPD(bnet, index.do(2),    'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.u(2)}     = gaussian_CPD(bnet, index.u(2),     'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));
bnet.CPD{index.f(2)}     = gaussian_CPD(bnet, index.f(2),     'mean', zeros(3,1),    'cov', sUnknown.*eye(d,d));

