
function bnet = buildGaussSum(W, mu, Sigma)

N = 3;
x1 = 1;  x2 = 2;  x3 = 3;  
dag = zeros(N,N);
dag(1,3) = 1;
dag(2,3) = 1;

ns = [1 1 1]; % no vector-valued

i_obs    = [x1 x3];
i_hidden = setdiff(1:N, i_obs);

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes, 'observed', i_obs);

s1 = Sigma(1);  m1  = mu(1);
s2 = Sigma(2);  m2  = mu(2);
s3 = Sigma(3);  m3  = mu(3);

m = [m1 m2 m3];
s = [s1 s2 s3];

for i=1:N
    if i == 3
        % x3 = x1 + x2
        bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', m(i), 'cov', s(i), 'weights', W, 'clamp_mean', 1, 'clamp_weights', 1);
    else
        bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', m(i), 'cov', s(i), 'clamp_mean', 1, 'clamp_weights', 1);
    end
end