clear
close all
clc

max_NB  = 20;
a_PLUQABA   = zeros(max_NB,1);
m_PLUQABA   = zeros(max_NB,1);
t_PLUQABA   = zeros(max_NB,1);
a_LUABA     = zeros(max_NB,1);
m_LUABA     = zeros(max_NB,1);
t_LUABA     = zeros(max_NB,1);

for NB = 4 : max_NB
    S_dmodel  = 1e-2;
    S_ymodel  = 1e-4;
    
    % Build the articulated chain
    dmodel   = autoTree(NB, 1, 1.5, 1);
    dmodel   = autoTreeStochastic(dmodel);
    dmodel.gravity = [0; -9.81; 0];
    
    % Build the MAP problem to solve the forward dynamics
    ymodel_LU  = autoSensSIEfriction(dmodel);
    ymodel_LU  = autoSensStochastic( ymodel_LU );
    mySens     = sensors( ymodel_LU );
    dmodel_LU  = autoTreeStochastic(dmodel);
    % set Sw_inv to zero
    i = dmodel_LU.Sw_inv.i;
    j = dmodel_LU.Sw_inv.j;
    for k = 1:length(i)
        dmodel_LU.Sw_inv = set(dmodel_LU.Sw_inv, zeros(size(dmodel_LU.Sw_inv(i(k),j(k)))), i(k), j(k));
    end
    % set Sv_inv to zero
    i = dmodel_LU.Sv_inv.i;
    j = dmodel_LU.Sv_inv.j;
    for k = 1:length(i)
        dmodel_LU.Sv_inv = set(dmodel_LU.Sv_inv, eye(size(dmodel_LU.Sv_inv(i(k),j(k)))), i(k), j(k));
    end
    
    %%
    %  Numerically check if there exist P and Q in the LU factorisation:
    %
    %             P*[D(q);Y(q)]*Q = L*U
    %
    %  which do not depend on thespecific q.
    
    myModel   = model(dmodel_LU);
    myLU      = LUABA(myModel, mySens);
    
    q     = rand(dmodel.NB,1)*10;
    dq    = rand(dmodel.NB,1);
    myLU  = myLU.setState(q,dq);
    D     = sparse(myLU.iDs, myLU.jDs, myLU.Ds, 19*dmodel.NB, 26*dmodel.NB);
    Y     = myLU.IDsens.sensorsParams.Ys;
    b     = [sparse(myLU.ibs, ones(size(myLU.ibs)), myLU.bs, 19*NB, 1); ones(7*NB,1)];

    % Apply the LU factorisation without the permutaations
    % that lead to the minimum number of fill-in.
    [L,U,P] = lu([D; Y]);
    
    [a1_min,m1_min] = lu_cost(P*[D; Y]);
    [a1_fb ,m1_fb]  = fb_cost(L, U, b);
    a_LUABA(NB)     = a1_min + a1_fb;
    m_LUABA(NB)     = m1_min + m1_fb;
    
    A = P*[D; Y];
    tstart = tic;
    [L,U] = lu(A);
    L\(U\b); 
    t_LUABA(NB) = toc(tstart);
    
    % Apply the LU factorisation with the permutaations
    % that lead to the minimum number of fill-in.
    [~,~,P,Q] = lu([D; Y]);
    
    [minL, minU]  = lu(P*[D; Y]*Q, 0);
    [a1_min,m1_min] = lu_cost(P*[D; Y]*Q);
    [a1_fb ,m1_fb]  = fb_cost(minL, minU, b);
    a_PLUQABA(NB)     = a1_min + a1_fb;
    m_PLUQABA(NB)     = m1_min + m1_fb;
    
    A = P*[D; Y]*Q;
    tstart = tic;
    [L,U] = lu(A);
    L\(U\b); 
    t_PLUQABA(NB) = toc(tstart);

end

fH = figure;
fontSize = 12;
axes('Units', 'normalized', 'Parent',fH, 'FontSize', fontSize);
hold on;
plot(8 : max_NB, [a_LUABA(8 : max_NB) a_PLUQABA(8 : max_NB)], 'LineWidth', 6)
grid on;
legend({'Forward dynamics solved with LU factorization', 'Forward dynamics solved with permuted LU factorization'}, 'Interpreter', 'latex', 'FontSize', fontSize,  'Location','northWest','Orientation','vertical');
ylabel({'Number of floating point additions'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'DefaultAxesFontName', 'Times New Roman');
xlabel('$N_B$ (number of links)', 'Interpreter', 'latex', 'FontSize', fontSize);

save2pdf('a_luSIE.pdf',fH,600);

fH = figure;
fontSize = 12;
axes('Units', 'normalized', 'Parent',fH, 'FontSize', fontSize);
hold on;
plot(8 : max_NB, [m_LUABA(8 : max_NB) m_PLUQABA(8 : max_NB)], 'LineWidth', 6)
grid on;
legend({'Forward dynamics solved with LU factorization', 'Forward dynamics solved with permuted LU factorization'}, 'Interpreter', 'latex', 'FontSize', fontSize,  'Location','northWest','Orientation','vertical');
ylabel({'Number of floating point multiplications'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'DefaultAxesFontName', 'Times New Roman');
xlabel('$N_B$ (number of links)', 'Interpreter', 'latex', 'FontSize', fontSize);

save2pdf('m_luSIE.pdf',fH,600);

fH = figure;
fontSize = 12;
axes('Units', 'normalized', 'Parent',fH, 'FontSize', fontSize);
hold on;
plot(8 : max_NB, [t_LUABA(8 : max_NB) t_PLUQABA(8 : max_NB)], 'LineWidth', 6)
grid on;
legend({'Forward dynamics solved with LU factorization', 'Forward dynamics solved with permuted LU factorization'}, 'Interpreter', 'latex', 'FontSize', fontSize,  'Location','northWest','Orientation','vertical');
ylabel({'Execution time [sec]'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'DefaultAxesFontName', 'Times New Roman');
xlabel('$N_B$ (number of links)', 'Interpreter', 'latex', 'FontSize', fontSize);

save2pdf('t_luSIE.pdf',fH,600);
