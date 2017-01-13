clear
close all
clc

max_NB  = 5;
max_y   = 20;
a_ABA   = zeros(max_NB,1);
m_ABA   = zeros(max_NB,1);
c_ABA   = zeros(max_NB,1);
a_map   = zeros(max_NB,max_y);
m_map   = zeros(max_NB,max_y);
c_map   = zeros(max_NB,max_y);
y_label = {'tau1'; 'd2q1'; 'f1'; 'fx1'; 'a1'; 'fB1'};
y_label = [y_label; 'tau2'; 'd2q2'; 'f2'; 'fx2'; 'a2'; 'fB2'];

for NB = 2 : max_NB
   % Build the articulated chain
   dmodel   = autoTree(NB, 1, 1, 1);
   dmodel   = autoTreeStochastic(dmodel);
   dmodel.gravity = [0; -9.81; 0];
   
   % Build the MAP problem to solve the forward dynamics
   ymodel_LU  = autoSensABA(dmodel);
   ymodel     = ymodel_LU;
   
   y_label = [y_label; ['tau' num2str(NB)]; ['d2q' num2str(NB)]; ['f' num2str(NB)]; ['fx' num2str(NB)]; ['a' num2str(NB)]; ['fB' num2str(NB)]];

   for ny = 1 : max_y
      r = randi([1 6*NB],1,ny);
      for k = 1 : ny
         ymodel     = addSens(ymodel, y_label{r(k)});
      end
      
      ymodel     = autoSensStochastic( ymodel );
      mySens     = sensors( ymodel );
      dmodel     = autoTreeStochastic(dmodel);
      % set Sw_inv to zero
      i = dmodel.Sw_inv.i;
      j = dmodel.Sw_inv.j;
      for k = 1:length(i)
         dmodel.Sw_inv = set(dmodel.Sw_inv, zeros(size(dmodel.Sw_inv(i(k),j(k)))), i(k), j(k));
      end
      % set Sv_inv to zero
      i = dmodel.Sv_inv.i;
      j = dmodel.Sv_inv.j;
      for k = 1:length(i)
         dmodel.Sv_inv = set(dmodel.Sv_inv, eye(size(dmodel.Sv_inv(i(k),j(k)))), i(k), j(k));
      end
      
      
      myModel = model(dmodel);
      myMAP   = MAP(myModel, mySens);
      q       = rand(dmodel.NB,1)*10;
      dq      = rand(dmodel.NB,1);
      y       = rand(ymodel.m,1);
      myMAP   = myMAP.setState(q,dq);
      myMAP   = myMAP.setY(y);
      myMAP   = myMAP.solveID();
      
      D = sparse(myMAP.iDs, myMAP.jDs, myMAP.Ds, 19*NB, 26*NB);
      Y = myMAP.IDsens.sensorsParams.Ys;
      b = [-sparse(myMAP.ibs, ones(size(myMAP.ibs)), myMAP.bs, 19*NB, 1); y];
      d = (D'*D + Y'*Y)\([D' Y']*b);
      d = d(myMAP.id,1);
      
      if norm(d-myMAP.d) > 1e-12
         error(['The MAP computations are not consistent with the least-squares solution. Error is: ', num2str(norm(d-myMAP.d))])
      end
      
      if ny == 1
         Y_ABA = Y(1:end-ymodel.sizes{end}, :);
         
         a_ABA(NB,1) = 205*NB - 248;
         m_ABA(NB,1) = 224*NB - 259;
         c_ABA(NB,1) = condest([D; Y_ABA]);
      end
      
      % Apply the LU facorisation and compute its cost.
      % Pvt is the pivot matrix which is shown to coincide
      % with the mass matrix
      [U,p,S] = chol(D'*D + Y'*Y);
      if p ~= 0
         error('Unable to perform the Cholesky factorisation');
      end
      
      L = my_chol(S'*(D'*D + Y'*Y)*S);
      [ atm, mtm ] = mult_cost( [D; Y]', [D; Y]);
      [ ach, mch ] = chol_cost( S'*(D'*D + Y'*Y)*S );
      [ afb, mfb ] = fb_cost( L, L', [D' Y']*b );
      
      
      a_map(NB,ny) = ach + afb + atm;
      m_map(NB,ny) = mch + mfb + mtm;
      c_map(NB,ny) = sqrt(condest([D; Y]'*[D; Y]));
   end
end

i_plot = round(linspace(1, max_y, 6));

plot([a_ABA, a_map(:,i_plot(1)), a_map(:,i_plot(2)), a_map(:,i_plot(3)), a_map(:,i_plot(4)), a_map(:,i_plot(5)), a_map(:,i_plot(6)) ])
legend('ABA', ['chol MAP' num2str(i_plot(1))],  ['chol MAP' num2str(i_plot(2))],  ['chol MAP' num2str(i_plot(3))], ['chol MAP' num2str(i_plot(4))],  ['chol MAP' num2str(i_plot(5))], ['chol MAP' num2str(i_plot(6))])

figure

plot([c_ABA, c_map(:,i_plot(1)), c_map(:,i_plot(2)), c_map(:,i_plot(3)), c_map(:,i_plot(4)), c_map(:,i_plot(5)), c_map(:,i_plot(6)) ])
legend('ABA', ['cond MAP' num2str(i_plot(1))],  ['cond MAP' num2str(i_plot(2))],  ['cond MAP' num2str(i_plot(3))], ['cond MAP' num2str(i_plot(4))],  ['cond MAP' num2str(i_plot(5))], ['cond MAP' num2str(i_plot(6))])
