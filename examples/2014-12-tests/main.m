clc
close all
clear all

global Sm Su
global dSv_inv idSv_inv jdSv_inv
global dSw_inv idSw_inv jdSw_inv
global dSy_inv idSy_inv jdSy_inv
sModel  = 1e1;
sMeas   = 1e-5;
sUknown = 1e3;


nbVals = logspace(log10(1), log10(50), 10);
nbVals = round(nbVals);
nbVals = unique(nbVals','rows')';

for nb = nbVals
   NB = nb;
   model = autoTree(NB, 1+rand(1)*ceil(NB/4));
   NB = NB + 5;
   n  = NB;
   idSw_inv = []; jdSw_inv = []; dSw_inv=[];
   idSv_inv = []; jdSv_inv = []; dSv_inv=[];
   for i = 1 : NB
      Sm.a{i}   = sModel.*generateSPDmatrix(6);
      Sm.fB{i}  = sModel.*generateSPDmatrix(6);
      Sm.f{i}   = sModel.*generateSPDmatrix(6);
      Sm.tau{i} = sModel.*generateSPDmatrix(1);
      
      [ii, jj, ss] = submatrixSparse(19*i-18, 19*i-18, ...
         [inv(Sm.a{i})  zeros(6,6)    zeros(6,6)   zeros(6,1); ...
         zeros(6,6)     inv(Sm.fB{i}) zeros(6,6)   zeros(6,1); ...
         zeros(6,6)     zeros(6,6)    inv(Sm.f{i}) zeros(6,1); ...
         zeros(1,6)     zeros(1,6)    zeros(1,6)   inv(Sm.tau{i})]);
      idSv_inv = [idSv_inv; ii];
      jdSv_inv = [jdSv_inv; jj];
      dSv_inv  = [dSv_inv;  ss];
       
      Su.fx{i}  = sUknown.*generateSPDmatrix(6);
      Su.d2q{i} = sUknown.*generateSPDmatrix(1);
     
      %dSv_inv((i-1)*19+(1:19), 1)   = [...
      %   1./diag(Sm.a{i});  1./diag(Sm.fB{i}); ... 
      %   1./diag(Sm.f{i});  1./diag(Sm.tau{i})];
      
      [ii, jj, ss] = submatrixSparse(7*i-6, 7*i-6, [inv(Su.fx{i}) zeros(6,1); zeros(1,6) inv(Su.d2q{i})]);
      idSw_inv = [idSw_inv; ii];
      jdSw_inv = [jdSw_inv; jj];
      dSw_inv  = [dSw_inv;  ss];
      %dSw_inv((i-1)*7+(1:7), 1)   = [...
      %   1./diag(Su.fx{i});  1./diag(Su.d2q{i})];      
   end   
   model = floatbase(model);
   fx = cell(NB);
   for i = 1 : NB
      fx{i} = rand(6,1);
   end
   q   = rand(NB,1);
   dq  = rand(NB,1);
   d2q = rand(NB,1);
   
   % build the measurements
   
   ny = 0;
   my = 1;
   ry = randfixedsum(NB, 1, 4*NB, 2, 6);
   ry = round(ry);
   idSy_inv = []; jdSy_inv = []; dSy_inv=[];
   for i = 1 : NB
      if ry(i) >= 6
         ny = ny + 1;
         y.sizes{ny,1}  = 6;
         y.labels{ny,1} = ['a' num2str(i)];
         y.S{ny,1}      = sMeas.*eye(6);
         [ii, jj, ss] = submatrixSparse(my, my, inv(y.S{ny,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 6;
      end
      if ry(i) >= 5
         ny = ny + 1;
         y.sizes{ny,1} = 6;
         y.labels{ny,1} = ['fB' num2str(i)];
         y.S{ny,1}      = sMeas.*eye(6);
         [ii, jj, ss] = submatrixSparse(my, my, inv(y.S{ny,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 6;         
      end
      if ry(i) >= 4
         ny = ny + 1;
         y.sizes{ny,1} = 6;
         y.labels{ny,1} = ['f' num2str(i)];
         y.S{ny,1}      = sMeas.*eye(6);
         [ii, jj, ss] = submatrixSparse(my, my, inv(y.S{ny,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 6;         
      end
      if ry(i) >= 3
         ny = ny + 1;
         y.sizes{ny,1} = 1;
         y.labels{ny,1} = ['tau' num2str(i)];
         y.S{ny,1}      = sMeas.*eye(1);
         [ii, jj, ss] = submatrixSparse(my, my, inv(y.S{ny,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 1;         
      end
      if ry(i) >= 2
         ny = ny + 1;
         y.sizes{ny,1} = 6;
         y.labels{ny,1} = ['fx'  num2str(i)];
         y.S{ny,1}      = sMeas.*eye(6);
         [ii, jj, ss] = submatrixSparse(my, my, inv(y.S{ny,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 6;         
      end
      if ry(i) >= 1
         ny = ny + 1;
         y.sizes{ny,1} = 1;
         y.labels{ny,1} = ['d2q' num2str(i)];
         y.S{ny,1}      = sMeas.*eye(1);
         [ii, jj, ss] = submatrixSparse(my, my, inv(y.S{ny,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 1;         
      end
   end
   
   % solve ID with the RNEA
   [tau_rne, a_rne, fB_rne, f_rne]        =  ID(model, q, dq, d2q, fx);
   
   % solve redundant ID with matrix inversion
   [y, Y]          = calcMeas (y, a_rne, fB_rne, f_rne, d2q, fx, ny, tau_rne, NB);
   [d_mat, Sd_mat] = mID(model, q, dq, Y, y);
   
   % solve redundant ID with Cholesky
   [ys, Ys, ym, Yx, Yy] = cholMeas (y, ny, NB);   
   [d_cho, Sd_cho, S]   = cID(model, q, dq, Ys, ys);
   
   % solve redundant ID without matrix inversion
   sparseModel = calcSparse(model);
   d_sol       = sID(model, q, dq, Yx, Yy, ym, S, sparseModel);
   
   % solve ID with a Bayesian network
   [bnet , y] = buildBnet(model, y);
   [d_net, Sd_net, bnet]  = nID(model, q, dq, d2q, fx, y, bnet);
   
   e_mat = max(abs(d_net - d_mat)) + max(diag(Sd_net - Sd_mat))
%    if (norm(e_mat) > 1e-5)
%       error('Something wrong with the mat method');
%    end

   e_sol = max(abs(d_net - d_sol))
%    if (norm(e_sol) > 1e-5)
%       error('Something wrong with the sol method');
%    end

   e_cho = max(abs(d_net - d_cho)) + max(diag(Sd_net - Sd_cho))
%    if (norm(e_cho) > 1e-5)
%       error('Something wrong with the chol method');
%    end
   
%    max(max(cell2mat(Sd_cho.Sf) - cell2mat(Sd_net.Sf))) + ...
%       max(max(cell2mat(Sd_cho.Sa) - cell2mat(Sd_net.Sa))) + ...
%       max(max(cell2mat(Sd_cho.Sfx) - cell2mat(Sd_net.Sfx))) + ...
%       max(max(cell2mat(Sd_cho.SfB) - cell2mat(Sd_net.SfB))) + ...
%       max(max(cell2mat(Sd_cho.Stau) - cell2mat(Sd_net.Stau))) + ...
%       max(max(cell2mat(Sd_cho.Sd2q) - cell2mat(Sd_net.Sd2q)))
end
