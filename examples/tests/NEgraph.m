clear all
close all
clc
%dim initialization
d = 3;   %vector dim
n = 2;   %number of links

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

k = n*6+4;

% indeces
i_p0 = 1;       i_do0 = 2;  
i_p  = [3 4];   i_do  = [5 6];
i_f  = [7 8 9]; i_u   = [10 11 12]; 
i_d2t= [13 14]; i_c= [15 16];

%uncertainties
sModel  = 2;
sUknown = 1;
sKnown  = 2;

% nodes dimesions
ns = [d*ones(1, 4*n+4), ones(1,n), d*ones(1, n)]; % all 3 dimensional vector for all but the d2t
dnodes = [];                      % no discrete nodes

% conditional relationship
dag = zeros(k);
dag(i_do0, i_do(1))   = 1;   dag(i_d2t(1), i_do(1)) = 1;
dag(i_do(1), i_do(2)) = 1;   dag(i_d2t(2), i_do(2)) = 1;

dag(i_p0, i_p(1))     = 1;   dag(i_do(1), i_p(1))   = 1;
dag(i_p(1), i_p(2))   = 1;   dag(i_do(2), i_p(2))   = 1;

dag(i_p(1), i_c(1))   = 1;   dag(i_do(1), i_c(1))   = 1;
dag(i_p(2), i_c(2))   = 1;   dag(i_do(2), i_c(2))   = 1;

dag(i_f(3), i_f(2))   = 1;   dag(i_c(2), i_f(2))    = 1;
dag(i_f(2), i_f(1))   = 1;   dag(i_c(1), i_f(1))    = 1;

dag(i_do(2), i_u(2))  = 1;   dag(i_f(2), i_u(2))    = 1;
dag(i_f(3), i_u(2))   = 1;   dag(i_u(3), i_u(2))    = 1;

dag(i_do(1), i_u(1))  = 1;   dag(i_f(1), i_u(1))    = 1;
dag(i_f(2), i_u(1))   = 1;   dag(i_u(2), i_u(1))    = 1;

%constant parameters
t1 = pi/4;    t2 = pi/4;
dt1 = 0.1;    dt2 =-0.1;
d2t1 = .1;    d2t2 =-.1;
l1 = 0.5;     l2 = 0.3;
lc1 = 0.25;   lc2 = 0.15;
g = 9.81;
m1 = 0.2;     m2 = 0.1;
I1z = 0.01;   I2z = 0.005;

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

%Initial conditions (forces and torques)

f3x=0; f3y=0; f3z=0;
u3x=0; u3y=0; u3z=0;

f3 = [f3x; f3y; f3z];
u3 = [u3x; u3y; u3z];

f{3} = f3;
u{3} = u3;

%Initial conditions (linear and angular accelerations)
o0x=0; o0y=0; o0z=0;
p0x=0; p0y=g; p0z=0;
do0x=0; do0y=0; do0z=10;

o0  = [o0x; o0y; o0z];
p0  = [p0x; p0y; p0z];
do0 = [do0x; do0y; do0z];

dt  = [dt1,  dt2];
d2t = [d2t1, d2t2];

% graph init
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

%graph topological reorder such that i-ancestors <  i-descendants
i_p0 = 1;       i_do0 = 2;  
i_p  = [3 4];   i_do  = [5 6];
i_f  = [7 8 9]; i_u   = [10 11 12]; 
i_d2t= [13 14]; i_c= [15 16];

i_p0    = find(bnet.order == i_p0);       i_do0 = find(bnet.order == i_do0);
i_p(1)  = find(bnet.order == i_p(1));    i_p(2) = find(bnet.order == i_p(2));
i_do(1) = find(bnet.order == i_do(1));  i_do(2) = find(bnet.order == i_do(2));
i_f(1)  = find(bnet.order == i_f(1));    i_f(2) = find(bnet.order == i_f(2));    i_f(3) = find(bnet.order == i_f(3));
i_u(1)  = find(bnet.order == i_u(1));    i_u(2) = find(bnet.order == i_u(2));    i_u(3) = find(bnet.order == i_u(3));
i_d2t(1)= find(bnet.order == i_d2t(1));i_d2t(2) = find(bnet.order == i_d2t(2));
i_c(1)  = find(bnet.order == i_c(1));    i_c(2) = find(bnet.order == i_c(2));

identity = eye(k);
new_dag = dag*identity(:,bnet.order);
new_dag = identity(bnet.order, :)*new_dag;
new_ns = ns(bnet.order);

bnet = mk_bnet(new_dag, new_ns, 'discrete', dnodes);

bnet.CPD{i_do0}    = gaussian_CPD(bnet, i_do0,    'mean', do0,    'cov', 1/100.*sUknown.*eye(d,d));
bnet.CPD{i_p0}     = gaussian_CPD(bnet, i_p0,     'mean', p0,     'cov', 10.*sUknown.*eye(d,d));
bnet.CPD{i_d2t(1)} = gaussian_CPD(bnet, i_d2t(1), 'mean', d2t(1), 'cov', sUknown);
bnet.CPD{i_d2t(2)} = gaussian_CPD(bnet, i_d2t(2), 'mean', d2t(2), 'cov', sUknown);
bnet.CPD{i_f(3)}   = gaussian_CPD(bnet, i_f(3),   'mean', f3,     'cov', sUknown.*eye(d,d));
bnet.CPD{i_u(3)}   = gaussian_CPD(bnet, i_u(3),   'mean', u3,     'cov', 1/10 .*sUknown.*eye(d,d));

%Useful definitions
c1 = cos(t1);  c2 = cos(t2);
s1 = sin(t1);  s2 = sin(t2);

R1 = [c1 -s1 0;
      s1  c1 0;
      0    0 1];

R2 = [c2 -s2 0;
      s2  c2 0;
      0    0 1];

R3 = [1 0 0;
      0 1 0;
      0 0 1];  

%Kinematic recursion
for i = 1:2
    if i==1
        Ri=R1;
    else
        Ri=R2;
    end
    if i==1
        o(:,i)  = Ri'*(o0+dt(i)*z0);
        
        %% do(:,i) = Ri'*(do_prev+d2t(i)*z0+dt(i)*cross(o_prev,z0));
        b = Ri'*dt(i)*cross(o0,z0);
        if (i_d2t(i) < i_do0)
            W = [Ri'*z0 Ri'];
        else
            W = [Ri' Ri'*z0];
        end
        bnet.CPD{i_do(i)}  = gaussian_CPD(bnet, i_do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
        b = cross(o(:,i),cross(o(:,i),r0(:,i)));
        if (i_p0 < i_do(i))
            W = [Ri' -vec_hat(r0(:,i))];
        else
            W = [-vec_hat(r0(:,i)) Ri'];
        end
        bnet.CPD{i_p(i)}  = gaussian_CPD(bnet, i_p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %%
        % c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));        
        b = cross(o(:,i),cross(o(:,i),rc(:,i)));
        if (i_p(i) < i_do(i))
            W = [eye(d) -vec_hat(rc(:,i))];
        else
            W = [-vec_hat(rc(:,i)) eye(d)];
        end        
        bnet.CPD{i_c(i)}  = gaussian_CPD(bnet, i_c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
    else
        o(:,i)  = Ri'*(o(:,i-1)+dt(i)*z0);
        
        %% do(:,i) = Ri'*(do_prev+d2t(i)*z0+dt(i)*cross(o_prev,z0));
        b = Ri'*dt(i)*cross(o(:,i-1),z0);
        if (i_d2t(i) < i_do0)
            W = [Ri'*z0 Ri'];
        else
            W = [Ri' Ri'*z0];
        end
        bnet.CPD{i_do(i)}  = gaussian_CPD(bnet, i_do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
        b = cross(o(:,i),cross(o(:,i),r0(:,i)));
        if (i_p(i-1) < i_do(i))
            W = [Ri' -vec_hat(r0(:,i))];
        else
            W = [-vec_hat(r0(:,i)) Ri'];
        end
        bnet.CPD{i_p(i)}  = gaussian_CPD(bnet, i_p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %%
        % c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));        
        b = cross(o(:,i),cross(o(:,i),rc(:,i)));
        if (i_p(i) < i_do(i))
            W = [eye(d) -vec_hat(rc(:,i))];
        else
            W = [-vec_hat(rc(:,i)) eye(d)];
        end   
        bnet.CPD{i_c(i)}  = gaussian_CPD(bnet, i_c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
    end
        
end

% Dynamic recursion
for i = 2:-1:1
    if i==2
        Ri=R3;
        Ii = I2;
    else
        Ri=R2;
        Ii = I1;
    end
    %% f(:,i)  = Ri'*f(:,i+1) + m(i,1).*c(:,i);
    b = zeros(d,1);
    if (i_f(i+1) < i_c(i))
        W = [Ri' eye(d).*m(i,1)];
    else
        W = [eye(d).*m(i,1) Ri];
    end
    bnet.CPD{i_f(i)}  = gaussian_CPD(bnet, i_f(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
    
    %% u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
    %    cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i));
    b = cross(o(:,i), Ii*o(:,i));
    
    Wt(:, :, 1) = Ii;
    Wt(:, :, 2) = vec_hat(rc(:,i)+r0(:,i));
    Wt(:, :, 3) = -vec_hat(rc(:,i))*Ri;
    Wt(:, :, 4) = Ri;
    [y,ind] = sort([i_do(i), i_f(i), i_f(i+1),  i_u(i+1)]);
    W = [Wt(:,:,ind(1)) Wt(:,:,ind(2)) Wt(:,:,ind(3)) Wt(:,:,ind(4))]; 
    
    bnet.CPD{i_u(i)}  = gaussian_CPD(bnet, i_u(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
end

engine{1} = jtree_inf_engine(bnet);

% set conditional nodes
evidence = cell(1,k);
%evidence{i_do0} = do0;
%evidence{i_do(2)} = [0; 0; 10.0];
%evidence{i_p0}  = p0;
%evidence{i_d2t(1)}  = d2t(1);
%evidence{i_d2t(2)}  = d2t(2);
%evidence{i_f(3)}  = f3;
%evidence{i_u(3)}  = u3;
[engine, ll] = enter_evidence(engine{1}, evidence);

marg = marginal_nodes(engine, [i_do(1)]);
doMu(:,1) = reshape(marg.mu,3,1);

% marg = marginal_nodes(engine, [i_do(2)]);
% doMu(:,2) = reshape(marg.mu,3,1);

marg = marginal_nodes(engine, [i_p(1)]);
pMu(:,1) = reshape(marg.mu,3,1);

marg = marginal_nodes(engine, [i_p(2)]);
pMu(:,2) = reshape(marg.mu,3,1);

marg = marginal_nodes(engine, [i_f(1)]);
fMu(:,1) = reshape(marg.mu,3,1);

marg = marginal_nodes(engine, [i_f(2)]);
fMu(:,2) = reshape(marg.mu,3,1);

marg = marginal_nodes(engine, [i_u(1)]);
uMu(:,1) = reshape(marg.mu,3,1);

marg = marginal_nodes(engine, [i_u(2)]);
uMu(:,2) = reshape(marg.mu,3,1);

marg = marginal_nodes(engine, [i_c(1)]);
cMu(:,1) = reshape(marg.mu,3,1);

marg = marginal_nodes(engine, [i_c(2)]);
cMu(:,2) = reshape(marg.mu,3,1);





