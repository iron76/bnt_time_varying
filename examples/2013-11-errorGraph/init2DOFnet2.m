function [bnet, index] = init2DOFnet2(n, k, d, l)

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

% indeces
i_p0   = 1;       i_do0 = 2;  
i_p    = [3 4];   i_do  = [5 6];
i_f    = [7 8 9]; i_u   = [10 11 12]; 
i_d2t  = [13 14]; i_c   = [15 16];
%error equations
i_edo  = [17 18]; i_ef  = [23 24];
i_ep   = [19 20]; i_eu  = [25 26];
i_ec   = [21 22];


% nodes dimesions
ns = [d*ones(1, 4*n+4), ones(1,n), d*ones(1, n), d*ones(1, l) ]; % all 3 dimensional vector for all but the d2t
dnodes = [];                      % no discrete nodes

% conditional relationship
dag = zeros(k);
dag(i_do(1), i_edo(1)) = 1;   dag(i_do0  , i_edo(1)) = 1;   dag(i_d2t(1), i_edo(1)) = 1;
dag(i_do(2), i_edo(2)) = 1;   dag(i_do(1), i_edo(2)) = 1;   dag(i_d2t(2), i_edo(2)) = 1;

dag(i_p(1), i_ep(1))  = 1;   dag(i_p0  , i_ep(1))    = 1;   dag(i_do(1), i_ep(1))   = 1;
dag(i_p(2), i_ep(2))  = 1;   dag(i_p(1), i_ep(2))    = 1;   dag(i_do(2), i_ep(2))   = 1;

dag(i_c(1), i_ec(1))   = 1;   dag(i_p(1), i_ec(1))   = 1;   dag(i_do(1), i_ec(1))   = 1;
dag(i_c(2), i_ec(2))   = 1;   dag(i_p(2), i_ec(2))   = 1;   dag(i_do(2), i_ec(2))   = 1;

dag(i_f(2), i_ef(2))   = 1;   dag(i_f(3), i_ef(2))   = 1;   dag(i_c(2), i_ef(2))    = 1;
dag(i_f(1), i_ef(1))   = 1;   dag(i_f(2), i_ef(1))   = 1;   dag(i_c(1), i_ef(1))    = 1;

dag(i_u(2), i_eu(2))  = 1;   dag(i_do(2), i_eu(2))   = 1;   dag(i_f(2), i_eu(2))    = 1;
dag(i_f(3), i_eu(2))  = 1;   dag(i_u(3) , i_eu(2))   = 1;

dag(i_u(1),  i_eu(1))  = 1;   dag(i_do(1), i_eu(1))  = 1;   dag(i_f(1), i_eu(1))    = 1;
dag(i_f(2),  i_eu(1))  = 1;   dag(i_u(2),  i_eu(1))  = 1;

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

i_edo(1) = find(bnet.order == i_edo(1));  i_edo(2) = find(bnet.order == i_edo(2));
i_ep(1)  = find(bnet.order == i_ep(1));    i_ep(2) = find(bnet.order == i_ep(2));
i_ec(1)  = find(bnet.order == i_ec(1));    i_ec(2) = find(bnet.order == i_ec(2));
i_ef(1)  = find(bnet.order == i_ef(1));    i_ef(2) = find(bnet.order == i_ef(2));
i_eu(1)  = find(bnet.order == i_eu(1));    i_eu(2) = find(bnet.order == i_eu(2));


identity = eye(k);
new_dag = dag*identity(:,bnet.order);
new_dag = identity(bnet.order, :)*new_dag;
new_ns = ns(bnet.order);

bnet = mk_bnet(new_dag, new_ns, 'discrete', dnodes);
index.p0    = i_p0;      index.p  = i_p;
index.do0   = i_do0;     index.do = i_do;
index.f     = i_f;       index.u  = i_u;
index.d2t   = i_d2t;     index.c  = i_c;

index.ep  = i_ep;        index.edo = i_edo;
index.ef  = i_ef;        index.eu  = i_eu;
index.ec  = i_ec;
