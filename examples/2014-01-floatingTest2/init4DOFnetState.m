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

i_o    = [1 2 3 4];
m_o    = [5 6 7 8];

% state nodes [q; dq]
i_x    = 9;

% sizes
s_o    = [1 1 1 1];       s_x    = 1;

% nodes dimesions
ns = [ s_o s_o s_x]; % all 3 dimensional vector for all but the d2q
dnodes = [];                      % no discrete nodes

% conditional relationship
dag = zeros(k);   

dag(i_o(1), m_o(1) ) = 1;   
dag(i_o(2), m_o(2) ) = 1;   dag(i_o(1), i_o(2)   ) = 1;
dag(i_o(3), m_o(3) ) = 1;   dag(i_o(2), i_o(3)   ) = 1;
dag(i_o(4), m_o(4) ) = 1;   dag(i_o(3), i_o(4)   ) = 1;


for i = i_o(3)
    dag(i_x, i)  = 1;
end

% graph init
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

%graph topological reorder such that i-ancestors <  i-descendants


i_o(1)  = find(bnet.order == i_o(1));   
i_o(2) = find(bnet.order == i_o(2));
i_o(3)  = find(bnet.order == i_o(3));    
i_o(4) = find(bnet.order == i_o(4));
m_o(1)  = find(bnet.order == m_o(1));    
m_o(2) = find(bnet.order == m_o(2));
m_o(3)  = find(bnet.order == m_o(3));    
m_o(4) = find(bnet.order == m_o(4));
i_x     = find(bnet.order == i_x);

identity = eye(k);
new_dag = dag*identity(:,bnet.order);
new_dag = identity(bnet.order, :)*new_dag;
new_ns = ns(bnet.order);

bnet = mk_bnet(new_dag, new_ns, 'discrete', dnodes);
index.o     = i_o;       index.mo  = m_o;
index.x     = i_x;

%init distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnet.CPD{index.mo(1)}  = gaussian_CPD(bnet, index.mo(1),   'mean', zeros(d,1), 'cov', So1, 'clamp_mean', 1, 'clamp_weights', 1,                                    'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mo(2)}  = gaussian_CPD(bnet, index.mo(2),   'mean', zeros(d,1), 'cov', So2, 'clamp_mean', 1, 'clamp_weights', 1,                                    'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mo(3)}  = gaussian_CPD(bnet, index.mo(3),   'mean', zeros(d,1), 'cov', So3, 'clamp_mean', 1, 'clamp_weights', 1,                                    'clamp_cov', 1, 'weights', eye(d));
bnet.CPD{index.mo(4)}  = gaussian_CPD(bnet, index.mo(4),   'mean', zeros(d,1), 'cov', So4, 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight,                'weights', eye(d));



