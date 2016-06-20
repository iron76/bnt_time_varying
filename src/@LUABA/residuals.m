function e = residuals(obj, d, pA, pa)


% Author: Francesco Nori
% Genova, June 2016
%
% Given an forward-dynamics problem formulated as [D; Y]*d + [bD; -y] = 0
% with d organized as follows:
%
%   d   = [dx_1, dx_2, ..., dx_obj.IDstate.n, dy_1, dy_2, ..., dy_obj.IDstate.n]
%
%   where:
%
%   dx_i = [a_i, fB_i, f_i, tau_i]
%
%   dy_i = [fx_i, d2q_i]
%
% computes permutation matrices P and Q such that P*[D;Y]*Q = L with L
% lower triangular. 

NB = obj.IDmodel.modelParams.NB;

pt(obj.id) = 1:length(obj.id);
d          = d(pt);

for i = 1 : length(obj.jd2q)
   itau_a(1+7*(i-1):7*i,1) = [obj.itau(i); obj.ia(1+6*(i-1):6*i)];
end

for i = 1 : length(obj.jd2q)
   jd2q_a(1+7*(i-1):7*i,1) = [obj.jd2q(i); obj.ja(1+6*(i-1):6*i)];
end

D  = sparse(obj.iDs, obj.jDs, obj.Ds, 19*NB, 26*NB);
Y  = obj.IDsens.sensorsParams.Ys;

% Let's build the matrix W such that (if lambda_i ~=0):
%
%   pi_a = f_i - Ii_a a_lambda_i
%
% and (if lambda_i == 0)
%
%   pi_A = f_i - Ii_A a_i
%
% In the case lambda_i = i-1, lambda_1 = 0 we have:
%
%   | p1_A |   | f1 |     | a1 |  
%   | p2_a |   | f2 |     | a2 |  
%   |  ..  | = | .. | + W | .. |
%   | pN_a |   | fN |     | aN | 
%
%              | -I1_A     0   ...     0    0 |
%              | -I2_a     0   ...     0    0 |
%        W   = |     0 -I3_a   ...     0    0 |
%              |                              |
%              |     0     0   ... -In_a    0 |
%
%
Ws = zeros(6*NB, 6*NB);
for i = 1 : NB
   %if obj.IDmodel.modelParams.parent(i) == 0
      Ws(1+6*(i-1):6*i, 1+6*(i-1):6*i) = -obj.IA{i};
   %else
   %   j = obj.IDmodel.modelParams.parent(i);
   %   Ws(1+6*(i-1):6*i, 1+6*(j-1):6*j) = -obj.Ia{i}*obj.Xup{i};
   %end
end
% Let's build the matrix Wr such that:
%
%   fi_b = fi_B - Ii ai
%
%   | f1_b |   | f1_B |      | a1 |  
%   | f2_b |   | f2_B |      | a2 |  
%   |  ..  | = |  ..  | + Wr | .. |
%   | fN_b |   | fN_B |      | aN | 
%
%              | -I1    0   ...   0 |
%              |   0  -I2   ...   0 |
%        Wr  = |   0    0   ...   0 |
%              |                    |
%              |   0    0   ... -IN |
%
%
Wr = zeros(6*NB, 6*NB);
for i = 1 : NB
   Wr(1+6*(i-1):6*i, 1+6*(i-1):6*i) = -obj.IDmodel.modelParams.I{i};
end

W = eye(26*NB, 26*NB);
W(obj.jF, obj.ja) = Ws;
W(obj.jfB, obj.ja) = Wr;
W = sparse(W);

% jp = obj.jF;

D = D * W^(-1);
Y = Y * W^(-1);

pA = cell2mat(pA');

%dbar = [d(obj.jtau); d(obj.jfx); pAa(end:-1:1,1); d(jd2q_a,1); d(obj.jfB)];
%dref = [d(obj.jtau); d(obj.jfx); d(obj.jF(end:-1:1, 1)); d(jd2q_a,1); d(obj.jfB)];
dbar = d;
dbar(obj.jF) = pA;
for i = 1 : NB
   dbar(obj.jfB(1+6*(i-1):6*i,1)) = dbar(obj.jfB(1+6*(i-1):6*i,1)) - obj.IDmodel.modelParams.I{i}*dbar(obj.ja(1+6*(i-1):6*i,1));
end

e   = pA  - d(obj.jF) - Ws * d(obj.ja);
e   = dbar - W*d;

I = [obj.jtau; obj.jfx; obj.jF(end:-1:1, 1); jd2q_a; obj.jfB];
e   = dbar(I) - W(I,I)*d(I);

return




