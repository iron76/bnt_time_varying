function [h, Sh, index, R, engine] = build4DOFbnetState(bnet, engine, index, x, P, y, o0, p0, do0)

load columnImplicit.mat
defineConstState;

bnet.CPD{index.x} = gaussian_CPD(bnet, index.x, 'mean', x, 'cov', P, 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);

q  = x(1:end/2);
dq = x(end/2+1:end);

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
         

%Useful quantities for implicit NE
param = [l1 l2 m1 m2 r1 r2 g I1z I2z];
Ci  = CY(x, param); 
bi  = bY(x, param);
dC  = dC_dx(x, y, param);
db  = dC * x;

indC = 1;
%Kinematic recursion
for i = 1:4
    
    j = (1+d*(i-1)):(d*i);
    j_prev = (1+d*(i-2)):(d*(i-1));
    
    Ri = R(:,:,i);
    if i==1
        p_prev  = p0;
        do_prev = do0;
    else
        index_do_prev = index.do(i-1);
        index_p_prev  =  index.p(i-1);
    end
    
    if i == 3 || i == 4
        
        %% o(:,i)  = Ri'*(o_prev+dq(i)*z0);
        W1 = Ci(indC:indC+2, colum.o(j_prev));
        W2 = dC(indC:indC+2, :);
        b = bi(indC:indC+2) - db(indC:indC+2);
        W = commuteOnOrder2(W1, W2, index.o(i-1), index.x);
        bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;

%         if i == 3
%             W1(3,:)
%             W2(3,:)*x - x(7)
%             b(3) 
%         else
%             W1(3,:) - [0 0 1]
%             W2(3,:)*x - x(8)
%             b(3)
%         end
        
        %% do(:,i) = Ri'*(do_prev+d2q(i)*z0+dq(i)*cross(o_prev,z0));
        b  = bi(indC:indC+2) - db(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.d2t(i));
        W2 = Ci(indC:indC+2, colum.do(j_prev));
        W3 = dC(indC:indC+2, :);
        W = commuteOnOrder3(W1, W2, W3, index.d2q(i), index_do_prev, index.x);
        bnet.CPD{index.do(i)}  = gaussian_CPD(bnet, index.do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;

        %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
        b  = bi(indC:indC+2) - db(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.p(j_prev));
        W2 = Ci(indC:indC+2, colum.do(j));
        W3 = dC(indC:indC+2, :);
        W = commuteOnOrder3(W1, W2, W3, index_p_prev, index.do(i), index.x);
        bnet.CPD{index.p(i)}  = gaussian_CPD(bnet, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;
        
        
        %% c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));
        b  = bi(indC:indC+2) - db(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.p(j)) ;
        W2 = Ci(indC:indC+2, colum.do(j));
        W3 = dC(indC:indC+2, :);
        W = commuteOnOrder3(W1, W2, W3, index.p(i), index.do(i), index.x);
        bnet.CPD{index.c(i)}  = gaussian_CPD(bnet, index.c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;

    else
        %% o(:,i)  = Ri'*o_prev;
        if i ~= 1
            W1 = Ci(indC:indC+2, colum.o(j_prev));
            W2 = dC(indC:indC+2, :);
            b = bi(indC:indC+2) - db(indC:indC+2);
            W = commuteOnOrder2(W1, W2, index.o(i-1), index.x);
            bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        else
            W = dC(indC:indC+2, :);
            b = bi(indC:indC+2) - db(indC:indC+2);
            bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        end
        indC = indC+3;

        
        %% do(:,i) = Ri'*do_prev;
        if i == 1
            b = Ci(indC:indC+2, colum.do0)*do_prev + bi(indC:indC+2) - db(indC:indC+2);
            W = dC(indC:indC+2, :);
            bnet.CPD{index.do(i)}  = gaussian_CPD(bnet, index.do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
            indC = indC+3;
        else
            b = bi(indC:indC+2) - db(indC:indC+2);
            W1 = Ci(indC:indC+2, colum.do(j_prev));
            W2 = dC(indC:indC+2, :);
            W = commuteOnOrder2(W1, W2, index.do(i-1), index.x);            
            bnet.CPD{index.do(i)}  = gaussian_CPD(bnet, index.do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
            indC = indC+3;            
        end

        
        %% p(:,i)  = Ri'*(p_prev+d2q(i)*z0)+2*dt(i)*cross(o(:,i), Ri'*z0)+cross(do(:,i),ri0_i(:,i))+cross(o(:,i),cross(o(:,i),ri0_i(:,i)));
        
        if i == 1
            b  = Ci(indC:indC+2, colum.p0)*p_prev + bi(indC:indC+2) - db(indC:indC+2);
            W1 = Ci(indC:indC+2, colum.d2t(i));
            W2 = Ci(indC:indC+2, colum.do(j));
            W3 = dC(indC:indC+2, :);
            W = commuteOnOrder3(W1, W2, W3, index.d2q(i), index.do(i), index.x);
            bnet.CPD{index.p(i)}  = gaussian_CPD(bnet, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
            indC = indC+3;
        else
            b = bi(indC:indC+2) - db(indC:indC+2);
            W1 = Ci(indC:indC+2, colum.p(j_prev));
            W2 = Ci(indC:indC+2, colum.d2t(i));
            W3 = Ci(indC:indC+2, colum.do(j));
            W4 = dC(indC:indC+2, :);
            W = commuteOnOrder4(W1, W2, W3, W4, index_p_prev, index.d2q(i), index.do(i), index.x);
            bnet.CPD{index.p(i)}  = gaussian_CPD(bnet, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
            indC = indC+3;
        end
        
        
                
        %% c(:,i)  = p(:,i) + cross(do(:,i), riC_i(:,i)) +cross(o(:,i),cross(o(:,i),riC_i(:,i)));
        b = bi(indC:indC+2) - db(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.p(j));
        W2 = Ci(indC:indC+2, colum.do(j));
        W3 = dC(indC:indC+2, :);
        W = commuteOnOrder3(W1, W2, W3, index.p(i), index.do(i), index.x);
        bnet.CPD{index.c(i)}  = gaussian_CPD(bnet, index.c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W); 
        indC = indC+3;

    end
end

j = (1+d*4):(d*5);

% f(5) = Ree * fee;
b = bi(indC:indC+2) - db(indC:indC+2);
W1 = Ci(indC:indC+2, colum.fee);
W2 = dC(indC:indC+2, :);
W = commuteOnOrder2(W1, W2, index.fee, index.x);
bnet.CPD{index.f(5)}  = gaussian_CPD(bnet, index.f(5),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
indC = indC+3;

% u(5) = Ree * uee;
b = bi(indC:indC+2) - db(indC:indC+2);
W1 = Ci(indC:indC+2, colum.uee);
W2 = dC(indC:indC+2, :);
W = commuteOnOrder2(W1, W2, index.uee, index.x);
bnet.CPD{index.u(5)}  = gaussian_CPD(bnet, index.u(5),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
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
        b  = bi(indC:indC+2) - db(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.f(j_post));
        W2 = Ci(indC:indC+2, colum.c(j));
        W3 = Ci(indC:indC+2, colum.fbb);
        W4 = dC(indC:indC+2, :);
        W = commuteOnOrder4(W1, W2, W3, W4, index.f(i+1), index.c(i), index.fbb, index.x);
        bnet.CPD{index.f(i)}  = gaussian_CPD(bnet, index.f(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;
        
        
        %% u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
        %    cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i)) + cross(-rc(:,i)-r0(:,i), Rbb*fbb) + Rbb*ubb;
        b = bi(indC:indC+2) - db(indC:indC+2);
        Wt{1} = Ci(indC:indC+2, colum.do(j));
        Wt{2} = Ci(indC:indC+2, colum.f(j));
        Wt{3} = Ci(indC:indC+2, colum.f(j_post));
        Wt{4} = Ci(indC:indC+2, colum.u(j_post));
        Wt{5} = Ci(indC:indC+2, colum.fbb);
        Wt{6} = Ci(indC:indC+2, colum.ubb);
        Wt{7} = dC(indC:indC+2, :);
        [~,ind] = sort([index.do(i), index.f(i), index.f(i+1),  index.u(i+1), index.fbb, index.ubb, index.x]);
        W = [Wt{ind(1)} Wt{ind(2)} Wt{ind(3)} Wt{ind(4)} Wt{ind(5)} Wt{ind(6)} Wt{ind(7)}];
        
        bnet.CPD{index.u(i)}  = gaussian_CPD(bnet, index.u(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);                
        indC = indC+3;

        
    else
        %% f(:,i)  = Ri*f(:,i+1) + m(i,1).*c(:,i);
        b =  bi(indC:indC+2) - db(indC:indC+2);
        W1 = Ci(indC:indC+2, colum.f(j_post));
        W2 = Ci(indC:indC+2, colum.c(j));
        W3 = dC(indC:indC+2, :);
        W = commuteOnOrder3(W1, W2, W3, index.f(i+1), index.c(i), index.x);
        bnet.CPD{index.f(i)}  = gaussian_CPD(bnet, index.f(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        indC = indC+3;
        
        %% u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
        %    cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i))
        b = bi(indC:indC+2) - db(indC:indC+2);
        Wt{1} = Ci(indC:indC+2, colum.do(j));
        Wt{2} = Ci(indC:indC+2, colum.f(j));
        Wt{3} = Ci(indC:indC+2, colum.f(j_post));
        Wt{4} = Ci(indC:indC+2, colum.u(j_post));
        Wt{5} = dC(indC:indC+2, :);
        [~,ind] = sort([index.do(i), index.f(i), index.f(i+1),  index.u(i+1) index.x]);
        W = [Wt{ind(1)} Wt{ind(2)} Wt{ind(3)} Wt{ind(4)} Wt{ind(5)}];
        
        bnet.CPD{index.u(i)}  = gaussian_CPD(bnet, index.u(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);                
        indC = indC+3;
    end
end

if isempty(engine)
    engine = jtree_inf_engine(bnet);
else
    engine = bnet_to_engine(bnet, engine);
end

% set conditional nodes
evidence = cell(1,k);

%evidence{index.md2q(4)}  = y(colum.d2t(4));
evidence{index.mf(2)}    = y(colum.f(4:6));
evidence{index.mu(2)}    = y(colum.u(4:6));
evidence{index.mf(5)}    = y(colum.f(13:15));
evidence{index.mu(5)}    = y(colum.u(13:15));
evidence{index.mo(4)}    = y(colum.o(10:12));
evidence{index.mp(4)}    = y(colum.p(10:12));
evidence{index.fbb}      = y(colum.fbb);
evidence{index.ubb}      = y(colum.ubb);
evidence{index.mx}       = [x(4) x(8)]';

[engine, ll] = enter_evidence(engine, evidence);


for i = 1 : n
    marg = marginal_nodes(engine, index.p(i));
    h{index.p(i)}  = marg.mu;
    Sh{index.p(i)} = marg.Sigma;
end

for i = 1 : n
    marg = marginal_nodes(engine, index.do(i));
    h{index.do(i)}  = marg.mu;
    Sh{index.do(i)} = marg.Sigma;
end

for i = 1 : n
    marg = marginal_nodes(engine, index.o(i));
    h{index.o(i)}  = marg.mu;
    Sh{index.o(i)} = marg.Sigma;
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

marg = marginal_nodes(engine, index.x);
h{index.x}  = marg.mu;
Sh{index.x} = marg.Sigma;

