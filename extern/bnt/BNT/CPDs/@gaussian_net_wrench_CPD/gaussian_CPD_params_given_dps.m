function [m, C, W] = gaussian_CPD_params_given_dps(CPD, domain, evidence)
% GAUSSIAN_CPD_PARAMS_GIVEN_EV_ON_DPS Extract parameters given evidence on all discrete parents
% function [m, C, W] = gaussian_CPD_params_given_ev_on_dps(CPD, domain, evidence)

ps = domain(1:end-1);
dps = ps(CPD.dps);
if isempty(dps)
  [m,C,W] = gaussian_net_wrench_get_canonical_parameters(CPD);
end
