function [h, Sh, index, R] = build4DOFbnetState(bnet, index, x, P, o0, p0, do0)

load columnImplicit.mat
defineConstState;

%bnet.CPD{index.x} = gaussian_CPD(bnet, index.x, 'mean', x, 'cov', P, 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
bnet.CPD{index.x} = gaussian_CPD(bnet, index.x, 'mean', 0, 'cov', 1e6, 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);

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
    end
    
    if i == 3 || i == 4
        
        %% o(:,i)  = Ri'*(o_prev+dq(i)*z0);
        % bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', b, 'cov', diag([1e2 1e2 1e-6]), 'weights', W);
        if i == 3
             
             W = 1;
             W = [W, Ri(3,3)];
             bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', [0 ]', 'cov', diag([1e-4 ]), 'weights', W);
        else
             %W = zeros(1,8);
             %W(1, 8) = 1;
             W = [Ri(3,3)];
             bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', x(8), 'cov', diag([1e-4 ]), 'weights', W);            
        end
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
        

    else
        %% o(:,i)  = Ri'*o_prev;
        if i == 1

            % bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', b, 'cov', diag([1e2 1e2 1e-6]));
            bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', [0 ], 'cov', diag([1e-3 ]));            
        else
            % bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', b, 'cov', diag([1e2 1e2 1e-6]), 'weights', W);
            bnet.CPD{index.o(i)}  = gaussian_CPD(bnet, index.o(i),  'mean', [0 ], 'cov', diag([1e-3 ]), 'weights', eye(1));
        end
    end
end


engine = jtree_inf_engine(bnet);

% set conditional nodes
evidence = cell(1,k);

evidence{index.mo(4)}    = [0]';%y(colum.o(10:12));

[engine, ll] = enter_evidence(engine, evidence);

for i = 1 : n
    marg = marginal_nodes(engine, index.o(i));
    h{index.o(i)}  = marg.mu;
    Sh{index.o(i)} = marg.Sigma;
end
marg = marginal_nodes(engine, index.x);
h{index.x}  = marg.mu;
Sh{index.x} = marg.Sigma;

