function [ ymodel ] = autoSensDNEA( dmodel , ymodel_ini, mask_q, mask_dq )
%AUTOSENSSNEA Generates a random sensor distribution articulated rigid body.
%   This function generates a structure that contains all measurements
%   needed to perform sparse inverse dynanmic compuations on the supplied
%   articulated rigid body (dmodel) assuming mutiple (and possibly
%   redudant) senors. The dynamic model is assumed to be in
%   featherstone-like format. The struct has the following fields:
%
%       ny - the number of available sensors
%
%       NB - the number of rigid links in the articulated rigid body
%
%   labels - a (ny,1) cell array containing the type of sensor. Possible type of
%            sensors are listed in the following where the subscript i
%            refers to the i-th link of the articulated rigid body.
%                     - a_i :  spatial accelration
%                     - fB_i:  net spatial force
%                     - f_i:   net sptaial force with parent
%                     - tau_i: joint torque
%                     - fx_i:  external spatial force
%                     - d2q_i: joint acceleration
%                     - dq_i:  joint velocity
%                     -  q_i:  joint position
%
%    sizes - a (ny,1) cell array containing the size of each sensor
%
%        m - the lenght of measurement vetor y (sum of all sizes)
%
%        Y - a (ny,NB) cell array describing the link between the sensors
%            available sensors and the vector of variables d
%            describing the dynamics of the articulated rigid body. In
%            particulat we have:
%                       Y{i,j} d_j = y_i
%            being y_i the i-th measurement and d_j the following vector of
%            dynamic variables associated to the articulated rigid body:
%                       d_j = [a_j, fB_j, f_j, tau_j, fx_j, d2q_j; q_j; dq_j]
%
%       Ys - a sparse representation of the whole-matrix Y
%
% Author: Francesco Nori
% Genova, Dec 2014

ymodel = ymodel_ini;

ymodel.Ys  = [ymodel.Ys, sparse([],[],[],ymodel.m, 2*dmodel.NB,0)];
for i = 1 : ymodel.ny
   ymodel.Y{i,dmodel.NB+1} = zeros(size(ymodel.Y{i,1},1), dmodel.NB);
   ymodel.Y{i,dmodel.NB+2} = zeros(size(ymodel.Y{i,1},1), dmodel.NB);
end
   
if nargin == 4
   for i = 1 : length(mask_q)
      if mask_q(i) == 1
         yi = zeros(1, 2*dmodel.NB);
         yi(1,i) = 1;
         ymodel.Ys =      [ ymodel.Ys;  [sparse(zeros(1, 26*dmodel.NB)) sparse(yi)]];         
         for j = 1 : dmodel.NB
            ymodel.Y{ymodel.ny+1,j} = zeros(1, 26);
         end
         ymodel.Y{ymodel.ny+1,dmodel.NB+1} = yi(1, 1:dmodel.NB);
         ymodel.Y{ymodel.ny+1,dmodel.NB+2} = yi(1, dmodel.NB+1:end);
         ymodel.m  = ymodel.m  + 1;
         ymodel.ny = ymodel.ny + 1;
         ymodel.sizes{ymodel.ny,1} = 1;
         ymodel.labels{ymodel.ny,1} = ['y_q' num2str(i)];         
      end
   end
   for i = 1 : length(mask_dq)
      if mask_dq(i) == 1
         yi = zeros(1, 2*dmodel.NB);
         yi(1, dmodel.NB + i) = 1;
         ymodel.Ys = [ ymodel.Ys;  sparse(zeros(1, 26*dmodel.NB)) sparse(yi)];
         for j = 1 : dmodel.NB
            ymodel.Y{ymodel.ny+1,j} = zeros(1, 26);
         end
         ymodel.Y{ymodel.ny+1,dmodel.NB+1} = yi(1, 1:dmodel.NB);
         ymodel.Y{ymodel.ny+1,dmodel.NB+2} = yi(1, dmodel.NB+1:end);
         ymodel.m  = ymodel.m  + 1;
         ymodel.ny = ymodel.ny + 1;
         ymodel.sizes{ymodel.ny,1} = 1;
         ymodel.labels{ymodel.ny,1} = ['y_dq' num2str(i)];
      end
   end
end

%% Angular velocity measurements
% sample 20% angular 
ry = ceil(dmodel.NB*0.2);
for i = 1 : ry
   % select a random body from 1 to NB
   ymodel.ny = ymodel.ny + 1;
   ymodel.sizes{ymodel.ny,1} = 3;
   ymodel.m  = ymodel.m  + ymodel.sizes{ymodel.ny,1};
   py = randi(dmodel.NB, 1);
   ymodel.labels{ymodel.ny,1} = ['y_omega' num2str(py)];
   
   for j = 1 : dmodel.NB
      ymodel.Y{ymodel.ny,j} = zeros(ymodel.sizes{ymodel.ny,1}, 26);
   end
   ymodel.Y{ymodel.ny,dmodel.NB+1} = zeros(ymodel.sizes{ymodel.ny,1}, dmodel.NB);
   ymodel.Y{ymodel.ny,dmodel.NB+2} = zeros(ymodel.sizes{ymodel.ny,1}, dmodel.NB); 
   ymodel.Ys = [ ymodel.Ys;  sparse(zeros(ymodel.sizes{ymodel.ny,1}, 28*dmodel.NB)) ];
end   

% ymodel.NB = dmodel.NB;
% ny = 0;
% 
% ry = randfixedsum(ymodel.NB, 1, 4*ymodel.NB, 2, 6);
% ry = round(ry);
% for i = 1 : ymodel.NB
%    if ry(i) >= 6
%       ny = ny + 1;
%       ymodel.sizes{ny,1} = 6;
%       ymodel.labels{ny,1} = ['y_a' num2str(i)];
%    end
%    if ry(i) >= 5
%       ny = ny + 1;
%       ymodel.sizes{ny,1} = 6;
%       ymodel.labels{ny,1} = ['y_fB' num2str(i)];
%    end
%    if ry(i) >= 4
%       ny = ny + 1;
%       ymodel.sizes{ny,1} = 6;
%       ymodel.labels{ny,1} = ['y_f' num2str(i)];
%    end
%    if ry(i) >= 3
%       ny = ny + 1;
%       ymodel.sizes{ny,1} = 1;
%       ymodel.labels{ny,1} = ['y_tau' num2str(i)];
%    end
%    if ry(i) >= 2
%       ny = ny + 1;
%       ymodel.sizes{ny,1} = 6;
%       ymodel.labels{ny,1} = ['y_fx'  num2str(i)];
%    end
%    if ry(i) >= 1
%       ny = ny + 1;
%       ymodel.sizes{ny,1} = 1;
%       ymodel.labels{ny,1} = ['y_d2q' num2str(i)];
%    end
% end
% 
% ymodel.m  = sum(cell2mat(ymodel.sizes));
% ymodel.ny = ny;
% ymodel.Y  = cell(ny, ymodel.NB);
% ymodel.Ys = cell(ny, ymodel.NB);
% 
% for i = 1 : ymodel.ny
%    for j = 1 : ymodel.NB
%       my = ymodel.sizes{i,1};
%       ymodel.Y{i,j}  = zeros(my, 26);
%       ymodel.Ys{i,j} = zeros(my, 26);
%       if strcmp(ymodel.labels{i,1}, ['y_a' num2str(j)])
%          d = rand(6,1);
%          ymodel.Y{i,j}  = [diag(d) zeros(6,6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
%          ymodel.Ys{i,j} = sparse(1:6,1:6,d, 6, 26);
%       end
%       if strcmp(ymodel.labels{i,1}, ['y_fB' num2str(j)])
%          d = rand(6,1);
%          ymodel.Y{i,j}  = [zeros(6,6) diag(d) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
%          ymodel.Ys{i,j} = sparse(1:6,7:12,d, 6, 26);
%       end
%       if strcmp(ymodel.labels{i,1}, ['y_f' num2str(j)])
%          d = rand(6,1);
%          ymodel.Y{i,j}  = [zeros(6,6) zeros(6,6) diag(d) zeros(6, 1) zeros(6,6) zeros(6, 1)];
%          ymodel.Ys{i,j} = sparse(1:6,13:18,d, 6, 26);
%       end
%       if strcmp(ymodel.labels{i,1}, ['y_tau' num2str(j)])
%          d = rand(1,1);
%          ymodel.Y{i,j}  = [zeros(1,6) zeros(1,6) zeros(1,6) diag(d) zeros(1,6) zeros(1, 1)];
%          ymodel.Ys{i,j} = sparse(1,19,d,1,26);
%       end
%       if strcmp(ymodel.labels{i,1}, ['y_fx' num2str(j)])
%          d = ones(6,1)*10;
%          ymodel.Y{i,j}  = [zeros(6,6) zeros(6,6) zeros(6, 6) zeros(6,1) diag(d) zeros(6, 1)];
%          ymodel.Ys{i,j} = sparse(1:6,20:25,d, 6, 26);
%       end
%       if strcmp(ymodel.labels{i,1}, ['y_d2q' num2str(j)])
%          d = rand(1,1);
%          ymodel.Y{i,j}  = [zeros(1,6) zeros(1,6) zeros(1,6) zeros(1, 1) zeros(1,6) diag(d)];
%          ymodel.Ys{i,j} = sparse(1,26,d,1,26);
%       end
%    end
% end
% 
% for i = 1 : ymodel.ny
%   for j = 1 : ymodel.NB
%     Yx{i,j}  = ymodel.Ys{i,j}(:, 1:19);
%     Yy{i,j}  = ymodel.Ys{i,j}(:, 20:end);
%   end
% end
% 
% ymodel.Ys = cell2mat([Yx Yy]);
% end
% 
% function [x,v] = randfixedsum(n,m,s,a,b)
% 
% % [x,v] = randfixedsum(n,m,s,a,b)
% %
% %   This generates an n by m array x, each of whose m columns
% % contains n random values lying in the interval [a,b], but
% % subject to the condition that their sum be equal to s.  The
% % scalar value s must accordingly satisfy n*a <= s <= n*b.  The
% % distribution of values is uniform in the sense that it has the
% % conditional probability distribution of a uniform distribution
% % over the whole n-cube, given that the sum of the x's is s.
% %
% %   The scalar v, if requested, returns with the total
% % n-1 dimensional volume (content) of the subset satisfying
% % this condition.  Consequently if v, considered as a function
% % of s and divided by sqrt(n), is integrated with respect to s
% % from s = a to s = b, the result would necessarily be the
% % n-dimensional volume of the whole cube, namely (b-a)^n.
% %
% %   This algorithm does no "rejecting" on the sets of x's it
% % obtains.  It is designed to generate only those that satisfy all
% % the above conditions and to do so with a uniform distribution.
% % It accomplishes this by decomposing the space of all possible x
% % sets (columns) into n-1 dimensional simplexes.  (Line segments,
% % triangles, and tetrahedra, are one-, two-, and three-dimensional
% % examples of simplexes, respectively.)  It makes use of three
% % different sets of 'rand' variables, one to locate values
% % uniformly within each type of simplex, another to randomly
% % select representatives of each different type of simplex in
% % proportion to their volume, and a third to perform random
% % permutations to provide an even distribution of simplex choices
% % among like types.  For example, with n equal to 3 and s set at,
% % say, 40% of the way from a towards b, there will be 2 different
% % types of simplex, in this case triangles, each with its own
% % area, and 6 different versions of each from permutations, for
% % a total of 12 triangles, and these all fit together to form a
% % particular planar non-regular hexagon in 3 dimensions, with v
% % returned set equal to the hexagon's area.
% %
% % Roger Stafford - Jan. 19, 2006
% 
% % Check the arguments.
% if (m~=round(m))|(n~=round(n))|(m<0)|(n<1)
%    error('n must be a whole number and m a non-negative integer.')
% elseif (s<n*a)|(s>n*b)|(a>=b)
%    error('Inequalities n*a <= s <= n*b and a < b must hold.')
% end
% 
% % Rescale to a unit cube: 0 <= x(i) <= 1
% s = (s-n*a)/(b-a);
% 
% % Construct the transition probability table, t.
% % t(i,j) will be utilized only in the region where j <= i + 1.
% k = max(min(floor(s),n-1),0); % Must have 0 <= k <= n-1
% s = max(min(s,k+1),k); % Must have k <= s <= k+1
% s1 = s - [k:-1:k-n+1]; % s1 & s2 will never be negative
% s2 = [k+n:-1:k+1] - s;
% w = zeros(n,n+1); w(1,2) = realmax; % Scale for full 'double' range
% t = zeros(n-1,n);
% tiny = 2^(-1074); % The smallest positive matlab 'double' no.
% for i = 2:n
%    tmp1 = w(i-1,2:i+1).*s1(1:i)/i;
%    tmp2 = w(i-1,1:i).*s2(n-i+1:n)/i;
%    w(i,2:i+1) = tmp1 + tmp2;
%    tmp3 = w(i,2:i+1) + tiny; % In case tmp1 & tmp2 are both 0,
%    tmp4 = (s2(n-i+1:n) > s1(1:i)); % then t is 0 on left & 1 on right
%    t(i-1,1:i) = (tmp2./tmp3).*tmp4 + (1-tmp1./tmp3).*(~tmp4);
% end
% 
% % Derive the polytope volume v from the appropriate
% % element in the bottom row of w.
% v = n^(3/2)*(w(n,k+2)/realmax)*(b-a)^(n-1);
% 
% % Now compute the matrix x.
% x = zeros(n,m);
% if m == 0, return, end % If m is zero, quit with x = []
% rt = rand(n-1,m); % For random selection of simplex type
% rs = rand(n-1,m); % For random location within a simplex
% s = repmat(s,1,m);
% j = repmat(k+1,1,m); % For indexing in the t table
% sm = zeros(1,m); pr = ones(1,m); % Start with sum zero & product 1
% for i = n-1:-1:1  % Work backwards in the t table
%    e = (rt(n-i,:)<=t(i,j)); % Use rt to choose a transition
%    sx = rs(n-i,:).^(1/i); % Use rs to compute next simplex coord.
%    sm = sm + (1-sx).*pr.*s/(i+1); % Update sum
%    pr = sx.*pr; % Update product
%    x(n-i,:) = sm + pr.*e; % Calculate x using simplex coords.
%    s = s - e; j = j - e; % Transition adjustment
% end
% x(n,:) = sm + pr.*s; % Compute the last x
% 
% % Randomly permute the order in the columns of x and rescale.
% rp = rand(n,m); % Use rp to carry out a matrix 'randperm'
% [ig,p] = sort(rp); % The values placed in ig are ignored
% x = (b-a)*x(p+repmat([0:n:n*(m-1)],n,1))+a; % Permute & rescale x
% end
% 
% function sparseModel = calcSparse(model)
% 
% NB  = model.NB;
% 
% for i = 1:model.NB
%    [ XJ, S{i} ] = jcalc( model.jtype{i}, 0);
%    [~, jn{i}] = size(S{i});
% end
% 
% sparseModel.nb = 0;
% sparseModel.nD = zeros(NB, 1);
% for i = 1 : NB
%    % b1  = Xup{i}*(-a_grav);
%    % OR
%    % b1 = crm(v{i})*vJ;
%    sparseModel.nb = sparseModel.nb + 6;
%    % b2 = crf(v{i})*model.I{i}*v{i};
%    sparseModel.nb = sparseModel.nb + 6;
%    
%    % Dii = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
%    sparseModel.nD(i) = sparseModel.nD(i) + 6 + jn{i}*6;
%    % Dii = [model.I{i} -eye(6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})];
%    sparseModel.nD(i) = sparseModel.nD(i) + 36 + 6;
%    % Dii = [zeros(6,6) eye(6) -eye(6) zeros(6, jn{i}) -inv(Xa{i}') zeros(6, jn{i})];
%    sparseModel.nD(i) = sparseModel.nD(i) + 6 + 6 + 36;
%    % Dii = [zeros(jn{i}, 6) zeros(jn{i}, 6) S{i}' -eye(jn{i}) zeros(jn{i}, 6) zeros(jn{i}, jn{i})];
%    sparseModel.nD(i) = sparseModel.nD(i) + 6*jn{i} + jn{i};
%    % Dij = [ Xup{i} zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i}) zeros(12+jn{i}, 24+2*jn{i})];
%    if model.parent(i) ~= 0
%       sparseModel.nD(i) = sparseModel.nD(i) + 36;
%    end
%    
%    ind_j  = find(model.parent == i);
%    for j = ind_j
%       % Dc{i,j} = [ zeros(12, 24+2*jn{i})
%       %     zeros(6,6) zeros(6,6) Xup{j}' zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})
%       %     zeros(jn{i}, 24+2*jn{i})];
%       sparseModel.nD(i) = sparseModel.nD(i) + 36;
%    end
% end
% 
% for i = 1 : NB
%    for j = 1 : NB
%       isparseModel.nD = [(j-1)*19 (i-1)*19];
%       
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(1:6),isparseModel.nD(2)+(1:6));
%       aa = aa';   bb = bb';   sparseModel.iD11(1:36,i,j) = aa(:);   sparseModel.jD11(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(1:6),isparseModel.nD(2)+(7:12));
%       aa = aa';   bb = bb';   sparseModel.iD12(1:36,i,j) = aa(:);   sparseModel.jD12(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(1:6),isparseModel.nD(2)+(13:18));
%       aa = aa';   bb = bb';   sparseModel.iD13(1:36,i,j) = aa(:);   sparseModel.jD13(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(1:6),isparseModel.nD(2)+19);
%       aa = aa';   bb = bb';   sparseModel.iD14(1: 6,i,j) = aa(:);   sparseModel.jD14(1: 6,i,j) = bb(:);
%       
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(7:12),isparseModel.nD(2)+(1:6));
%       aa = aa';   bb = bb';   sparseModel.iD21(1:36,i,j) = aa(:);   sparseModel.jD21(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(7:12),isparseModel.nD(2)+(7:12));
%       aa = aa';   bb = bb';   sparseModel.iD22(1:36,i,j) = aa(:);   sparseModel.jD22(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(7:12),isparseModel.nD(2)+(13:18));
%       aa = aa';   bb = bb';   sparseModel.iD23(1:36,i,j) = aa(:);   sparseModel.jD23(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(7:12),isparseModel.nD(2)+19);
%       aa = aa';   bb = bb';   sparseModel.iD24(1: 6,i,j) = aa(:);   sparseModel.jD24(1: 6,i,j) = bb(:);
%       
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(13:18),isparseModel.nD(2)+(1:6));
%       aa = aa';   bb = bb';   sparseModel.iD31(1:36,i,j) = aa(:);   sparseModel.jD31(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(13:18),isparseModel.nD(2)+(7:12));
%       aa = aa';   bb = bb';   sparseModel.iD32(1:36,i,j) = aa(:);   sparseModel.jD32(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(13:18),isparseModel.nD(2)+(13:18));
%       aa = aa';   bb = bb';   sparseModel.iD33(1:36,i,j) = aa(:);   sparseModel.jD33(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(13:18),isparseModel.nD(2)+19);
%       aa = aa';   bb = bb';   sparseModel.iD34(1: 6,i,j) = aa(:);   sparseModel.jD34(1: 6,i,j) = bb(:);
%       
%       [aa, bb] = meshgrid(isparseModel.nD(1)+19,isparseModel.nD(2)+(1:6));
%       aa = aa';   bb = bb';   sparseModel.iD41(1: 6,i,j) = aa(:);   sparseModel.jD41(1: 6,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+19,isparseModel.nD(2)+(7:12));
%       aa = aa';   bb = bb';   sparseModel.iD42(1: 6,i,j) = aa(:);   sparseModel.jD42(1: 6,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+19,isparseModel.nD(2)+(13:18));
%       aa = aa';   bb = bb';   sparseModel.iD43(1: 6,i,j) = aa(:);   sparseModel.jD43(1: 6,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+19,isparseModel.nD(2)+19);
%       aa = aa';   bb = bb';   sparseModel.iD44(1: 1,i,j) = aa;      sparseModel.jD44(1: 1,i,j) = bb(:);
%       
%       isparseModel.nD = [(j-1)*19 19*NB+(i-1)*7];
%       
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(1:6),isparseModel.nD(2)+(1:6));
%       aa = aa';   bb = bb';   sparseModel.iD15(1:36,i,j) = aa(:);   sparseModel.jD15(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(1:6),isparseModel.nD(2)+7);
%       aa = aa';   bb = bb';   sparseModel.iD16(1: 6,i,j) = aa(:);   sparseModel.jD16(1: 6,i,j) = bb(:);
%       
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(7:12),isparseModel.nD(2)+(1:6));
%       aa = aa';   bb = bb';   sparseModel.iD25(1:36,i,j) = aa(:);   sparseModel.jD25(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(7:12),isparseModel.nD(2)+7);
%       aa = aa';   bb = bb';   sparseModel.iD26(1: 6,i,j) = aa(:);   sparseModel.jD26(1: 6,i,j) = bb(:);
%       
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(13:18),isparseModel.nD(2)+(1:6));
%       aa = aa';   bb = bb';   sparseModel.iD35(1:36,i,j) = aa(:);   sparseModel.jD35(1:36,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+(13:18),isparseModel.nD(2)+7);
%       aa = aa';   bb = bb';   sparseModel.iD36(1: 6,i,j) = aa(:);   sparseModel.jD36(1: 6,i,j) = bb(:);
%       
%       [aa, bb] = meshgrid(isparseModel.nD(1)+19,isparseModel.nD(2)+(1:6));
%       aa = aa';   bb = bb';   sparseModel.iD45(1: 6,i,j) = aa(:);   sparseModel.jD45(1: 6,i,j) = bb(:);
%       [aa, bb] = meshgrid(isparseModel.nD(1)+19,isparseModel.nD(2)+7);
%       aa = aa';   bb = bb';   sparseModel.iD46(1: 1,i,j) = aa(:);   sparseModel.jD46(1: 1,i,j) = bb(:);
%    end
% end
% 
% for i = 1 : NB
%    sparseModel.ind_j{i}  = find(model.parent == i);
% end
% 
% end
