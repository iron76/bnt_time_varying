close all
clear all

rand('seed',0);
NB = 100;
model = autoTree(NB, 2.5);
NB = NB + 5;
n  = NB;
model = floatbase(model);
fx = cell(NB);
for i = 1 : NB
  fx{i} = rand(6,1);
end
q   = rand(NB,1);
dq  = rand(NB,1);
d2q = rand(NB,1);
% 
% % build the Bayesian network
bnet = buildBnet(model);

% load bnet.mat
% solve ID with the RNEA
[tau_rne, a_rne, fB_rne, f_rne]        =  ID(model, q, dq, d2q, fx);

% solve redundant ID with matrix inversion
[tau_mat, a_mat, fB_mat, f_mat]        = mID(model, q, dq, d2q, fx);

% solve ID with a Bayesian network
[tau_net, a_net, fB_net, f_net, bnet]  = nID(model, q, dq, d2q, fx, bnet);

%%%% SECOND PASS %%%%%
for i = 1 : NB
  fx{i} = rand(6,1);
end
d2q = rand(NB,1);
dq  = rand(NB,1);
q   = rand(NB,1);

% solve ID with the RNEA
tic;
[tau_rne, a_rne, fB_rne, f_rne]        =  ID(model, q, dq, d2q, fx);
t_rne = toc;

% solve redundant ID with matrix inversion
tic;
[tau_mat, a_mat, fB_mat, f_mat]        = mID(model, q, dq, d2q, fx);
t_mat = toc;

% solve ID with a Bayesian network
tic;
[tau_net, a_net, fB_net, f_net, bnet]  = nID(model, q, dq, d2q, fx, bnet);
t_net = toc;

fprintf(1, 'total time is (rne, mat, net): %.4f[ms], %.4f[ms], %.4f[ms] \r',...
    t_rne*1000, t_mat*1000, t_net*1000);

% norm(cell2mat(a_rne)-a_mat')
% norm(cell2mat(fB_rne)-fB_mat')
% norm(cell2mat(f_rne)-f_mat')
% norm(tau_rne-tau_mat)
% 
% norm(cell2mat(a_rne)-a_net')
% norm(cell2mat(fB_rne)-fB_net')
% norm(cell2mat(f_rne)-f_net')
% norm(tau_rne-tau_net)
