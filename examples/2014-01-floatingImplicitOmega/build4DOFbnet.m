function [h, Sh, bnet, R, avgTime, engine] = build4DOFbnet(bnet, index, y, q, dq, o0, do0, p0, start, engine)


load columnImplicit.mat
defineConst;


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
    if i~=1
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
        W = commuteOnOrder2(W1, W2, index.d2t(i), index_do_prev);
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
            b = Ci(indC:indC+2, colum.do0) * do0;
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
            b  = bi(indC:indC+2) + Ci(indC:indC+2, colum.p0)*p0;
            W1 = Ci(indC:indC+2, colum.d2t(i));
            W2 = Ci(indC:indC+2, colum.do(j));
            W = commuteOnOrder2(W1, W2, index.d2t(i), index.do(i));
            bnet.CPD{index.p(i)}  = gaussian_CPD(bnet, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        else
            b  = bi(indC:indC+2);
            W1 = Ci(indC:indC+2, colum.p(j_prev));
            W2 = Ci(indC:indC+2, colum.d2t(i));
            W3 = Ci(indC:indC+2, colum.do(j));
            W = commuteOnOrder3(W1, W2, W3, index_p_prev, index.d2t(i), index.do(i));
            bnet.CPD{index.p(i)}  = gaussian_CPD(bnet, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        end
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

tic;
if start == 1
    engine = gaussian_inf_engine(bnet);
end

% set conditional nodes
evidence = cell(k,1);

evidence{index.mp(4), 1}    = y(1:3);
evidence{index.mo(4), 1}    = y(4:6);
evidence{index.mf(2), 1}    = R(1:3,1:3,3)*y(7:9);
evidence{index.mu(2), 1}    = R(1:3,1:3,3)*y(10:12);
evidence{index.mf(5), 1}    = y(13:15);
evidence{index.mu(5), 1}    = y(16:18);
evidence{index.md2t(3), 1}  = y(19);
evidence{index.md2t(4), 1}  = y(20);
evidence{index.fbb, 1}      = y(21:23);
evidence{index.ubb, 1}      = y(24:26);

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
    marg = marginal_nodes(engine, index.d2t(i));
    h{index.d2t(i)}  = marg.mu;
    Sh{index.d2t(i)} = marg.Sigma;
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

avgTime = toc;

