function [bnet, LL, engine] = learn_params_em_modified(engine, evidence, max_iter, thresh)
% LEARN_PARAMS_EM Set the parameters of each adjustable node to their ML/MAP values using batch EM.
% [bnet, LLtrace, engine] = learn_params_em(engine, data, max_iter, thresh)
%
% data{i,l} is the value of node i in case l, or [] if hidden.
%   Suppose you have L training cases in an O*L array, D, where O is the num observed
%   scalar nodes, and N is the total num nodes.
%   Then you can create 'data' as follows, where onodes is the index of the observable nodes:
%      data = cell(N, L);
%      data(onodes,:) = num2cell(D);
%   Of course it is possible for different sets of nodes to be observed in each case.
%
% We return the modified bnet and engine.
% To see the learned parameters for node i, use the construct
%   s = struct(bnet.CPD{i}); % violate object privacy
% LLtrace is the learning curve: the vector of log-likelihood scores at each iteration.
%
% max_iter specifies the maximum number of iterations. Default: 10.
%
% thresh specifies the thresold for stopping EM. Default: 1e-3.
% We stop when |f(t) - f(t-1)| / avg < threshold,
% where avg = (|f(t)| + |f(t-1)|)/2 and f is log lik.  

if nargin < 3, max_iter = 10; end
if nargin < 4, thresh = 1e-3; end

verbose = 1;

loglik = 0;
previous_loglik = -inf;
converged = 0;
num_iter = 1;
LL = [];

while ~converged & (num_iter <= max_iter)
  [engine, loglik] = EM_step_modified(engine, evidence);
  if verbose, fprintf('EM iteration %d, ll = %8.4f\n', num_iter, loglik); end
  num_iter = num_iter + 1;
  converged = em_converged(loglik, previous_loglik, thresh);
  previous_loglik = loglik;
  LL = [LL loglik];
end
if verbose, fprintf('\n'); end

[mEng, nEng] = size(engine);
for k = 1 : nEng
    bnet{k} = bnet_from_engine(engine{k});
end

%%%%%%%%%

function [engine, loglik] = EM_step_modified(engine, cases)

[mEng, nEng] = size(engine);
loglik = 0;
for k = 1 : nEng
    bnet{k} = bnet_from_engine(engine{k}); % engine contains the old params that are used for the E step
    CPDs = bnet{k}.CPD; % these are the new params that get maximized
    num_CPDs = length(CPDs);
    adjustable = zeros(1,num_CPDs);
    for e=1:num_CPDs
        adjustable(e) = adjustable_CPD(CPDs{e});
    end
    adj = find(adjustable);
    n = length(bnet{k}.dag);
    
    for e=adj(:)'
        CPDs{e} = reset_ess(CPDs{e});
    end
    
    ncases = size(cases, 2);
    
    evidence = cases(:,k);
    [engine{k}, ll] = enter_evidence(engine{k}, evidence);
    loglik = loglik + ll;
    hidden_bitv = zeros(1,n);
    hidden_bitv(isemptycell(evidence))=1;
    for i=1:n
        e = bnet{k}.equiv_class(i);
        if adjustable(e)
            fmarg = marginal_family(engine{k}, i);
            if k == 1 
                CPDs{e} = update_ess_modified(CPDs{e}, CPDs{e}, fmarg, evidence, bnet{k}.node_sizes, bnet{k}.cnodes, hidden_bitv);
                CPDp{e} = CPDs{e};
            else
                CPDs{e} = update_ess_modified(CPDs{e}, CPDp{e}, fmarg, evidence, bnet{k}.node_sizes, bnet{k}.cnodes, hidden_bitv);
                CPDp{e} = CPDs{e};
            end 
        end
    end
end

for e=adj(:)'
    CPDs{e} = maximize_params_modified(CPDs{e});
end

for k = 1 : nEng
    engine{k} = update_engine_modified(engine{k}, CPDs, bnet{k}.CPD);
end


