function [bnet, index] = init4DOFnetState(n, k, d)

% graph initialization. graph contains the following nodes:
%
% - linear and angular link accelerations p, do including base [2*(n+1)]
%     p0, do0, p = [p1, p2]; do = [do1, do2],
% - forces and torques including end effector f, u [2*(n+1)]
%     f = [f1, f2, f3];  u = [u1, u2, u3];
% - joint acceletations d2q1, d2q2 [n]
%     d2q1, d2q2
% - center of mass acceleration [n]
%     c = [c1, c2]

defineConstState

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

% state nodes [q; dq]
i_x    = 65;

% sizes
s_p    = [3 3 3 3];        s_do  = [3 3 3 3];
s_d2q  = [1 1 1 1];
s_c    = [3 3 3 3];
s_f    = [3 3 3 3 3];      s_u   = [3 3 3 3 3];
s_fbb  = 3;                s_ubb = 3;
s_fee  = 3;                s_uee = 3;
s_o    = [3 3 3 3];       s_x    = 8;

% nodes dimesions
ns = [s_p s_do s_d2q s_c s_p s_do s_d2q s_c s_f s_u s_fbb s_ubb s_fee s_uee s_f s_u s_o s_o s_x]; % all 3 dimensional vector for all but the d2q
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

%state connections
for i = [i_p i_do i_d2q i_c]
    dag(i_x, i)  = 1;
end

for i = [i_f i_u]
    dag(i_x, i)  = 1;
end

for i = [i_f i_u]
    dag(i_x, i)  = 1;
end

for i = i_o
    dag(i_x, i)  = 1;
end

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

i_x     = find(bnet.order == i_x);

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
index.x     = i_x;

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



