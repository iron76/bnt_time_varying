
function bnet = buildGaussSum(W, mu, Sigma)

N = 2;
x1 = 1;  x2 = 2;
dag = zeros(N,N);
dag(1,2) = 1;

ns = [1 1]; % no vector-valued

i_obs    = [x2];
i_hidden = setdiff(1:N, i_obs);

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes, 'observed', i_obs);

s1 = Sigma(1);  m1  = mu(1);
s2 = Sigma(2);  m2  = mu(2);

m = [m1 m2];
s = [s1 s2];


% x2 = W*x1 + w
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', m(2), 'cov', s(2), 'weights', W, 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', 0.000001);
bnet.CPD{1} = gaussian_CPD(bnet, 1, 'mean', m(1), 'cov', s(1), 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);