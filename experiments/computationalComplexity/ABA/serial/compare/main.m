clear
close all
clc

max_NB  = 30;
a_matrixABA = zeros(max_NB,1);
m_matrixABA = zeros(max_NB,1);
t_matrixABA = zeros(max_NB,1);
a_ABA   = zeros(max_NB,1);
m_ABA   = zeros(max_NB,1);
a_myLUABA    = zeros(max_NB,1);
m_myLUABA    = zeros(max_NB,1);
a_PLUQABA   = zeros(max_NB,1);
m_PLUQABA   = zeros(max_NB,1);
t_PLUQABA   = zeros(max_NB,1);
a_LUABA     = zeros(max_NB,1);
m_LUABA     = zeros(max_NB,1);
t_LUABA     = zeros(max_NB,1);

for NB = 2 : max_NB
    S_dmodel  = 1e-2;
    S_ymodel  = 1e-4;
    
    % Build the articulated chain
    dmodel   = autoTree(NB, 1, 1.5, 0.9);
    dmodel   = autoTreeStochastic(dmodel);
    dmodel.gravity = [0; -9.81; 0];
    
    % Build the MAP problem to solve the forward dynamics
    ymodel_LU  = autoSensABA(dmodel);
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
    
    %%
    %  Checks the computational complexity of the ABA algorithm reformulating
    %  the ABA as the following factorisation:
    %
    %             P_ABA*[WL*D*WR; Y]*Q_ABA
    %
    %  which should result in a lower triangular matrix.
    a_ABA(NB,1) = 205*NB - 248;
    m_ABA(NB,1) = 224*NB - 259;
    [P_ABA, Q_ABA, WL, WR] = factorize(myLU);
    
    
    D     = sparse(myLU.iDs, myLU.jDs, myLU.Ds, 19*NB, 26*NB);
    Y     = myLU.IDsens.sensorsParams.Ys;
    b     = [sparse(myLU.ibs, ones(size(myLU.ibs)), myLU.bs, 19*NB, 1); ones(7*NB,1)];
    A_ABA = tril(P_ABA*[WL{4}*WL{3}*WL{2}*WL{1}*D*WR; Y]*Q_ABA);
    
    [a1_myABA, m1_myABA] = fw_cost(A_ABA, b);
    [a2_myABA, m2_myABA] = mult_cost( WL{2}, WL{1});
    [a3_myABA, m3_myABA] = mult_cost( WL{3}, WL{2}*WL{1});
    [a4_myABA, m4_myABA] = mult_cost( WL{4}, WL{3}*WL{2}*WL{1});
    [a5_myABA, m5_myABA] = mult_cost( WL{4}*WL{3}*WL{2}*WL{1}, D);
    
    a_matrixABA(NB,1) = a1_myABA + a2_myABA + a3_myABA + a4_myABA + a5_myABA;
    m_matrixABA(NB,1) = m1_myABA + m2_myABA + m3_myABA + m4_myABA + m5_myABA;
    
    L = P_ABA*[WL{4}*WL{3}*WL{2}*WL{1}*D*WR; Y]*Q_ABA;
    tstart = tic;
    WL{4}*WL{3}*WL{2}*WL{1}*D*WR;
    L\b; 
    t_matrixABA(NB) = toc(tstart);
    
    %%
    %  Performs a shuffle of the columns and rows so as to obtain a
    %  matrix on which we can gurantee the existence of the LU factorisation
    %  and proper defined pivots.
    
    B1 = D;
    B2 = B1(myLU.itau, :);
    B1(myLU.itau, : ) = [];
    B  = [B1; B2];
    
    A1 = [Y(:,[myLU.jfx; myLU.jtau]); B(:,[myLU.jfx; myLU.jtau])];
    A2 = [Y; B];
    A2(:, [myLU.jfx; myLU.jtau]) = [];
    A = [A1, A2];
    
    % Apply the LU facorisation and compute its cost.
    % Pvt is the pivot matrix which is shown to coincide
    % with the mass matrix
    [myL,myU,Pvt] = my_lu(A);
    [a1_my, m1_my]   = lu_cost(A);
    [a1_fb, m1_fb]  = fb_cost(myL, myU, b);
    a_myLUABA(NB,1)      = a1_my + a1_fb;
    m_myLUABA(NB,1)      = m1_my + m1_fb;
    
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

% plot([a_ABA, a_matrixABA, a_myLUABA, a_PLUQABA])
% legend('ABA', 'matrix ABA', 'LU factorisation ABA', 'P*LU*Q ABA')

fH = figure;
fontSize = 12;
axes('Units', 'normalized', 'Parent',fH, 'FontSize', fontSize);
hold on;
plot([a_LUABA, a_PLUQABA, a_matrixABA, a_ABA], 'LineWidth', 6)
grid on;
legend({'Forward dynamics solved with LU factorization', ...
    'Forward dynamics solved with permuted LU factorization', ...
    'Forward dynamics solved with forward substitution', ...
    'Forward dynamics solved with ABA'}, 'Interpreter', 'latex', 'FontSize', fontSize,  'Location','northWest','Orientation','vertical');
ylabel({'Number of floating point multiplications and additions'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'DefaultAxesFontName', 'Times New Roman');
xlabel('$N_B$ (number of links)', 'Interpreter', 'latex', 'FontSize', fontSize);

save2pdf('a_luAB.pdf',fH,600);

fH = figure;
fontSize = 12;
axes('Units', 'normalized', 'Parent',fH, 'FontSize', fontSize);
hold on;
plot([m_LUABA, m_PLUQABA, m_matrixABA, m_ABA], 'LineWidth', 6)
grid on;
legend({'Forward dynamics solved with LU factorization', ...
    'Forward dynamics solved with permuted LU factorization', ...
    'Forward dynamics solved with forward substitution', ...
    'Forward dynamics solved with ABA'}, 'Interpreter', 'latex', 'FontSize', fontSize,  'Location','northWest','Orientation','vertical');
ylabel({'Number of floating point multiplications and additions'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'DefaultAxesFontName', 'Times New Roman');
xlabel('$N_B$ (number of links)', 'Interpreter', 'latex', 'FontSize', fontSize);

save2pdf('m_luAB.pdf',fH,600);

fH = figure;
fontSize = 12;
axes('Units', 'normalized', 'Parent',fH, 'FontSize', fontSize);
hold on;
plot([m_LUABA+a_LUABA, m_PLUQABA+a_PLUQABA, m_matrixABA+a_matrixABA, m_ABA+a_ABA], 'LineWidth', 6)
grid on;
legend({'Forward dynamics solved with LU factorization', ...
    'Forward dynamics solved with permuted LU factorization', ...
    'Forward dynamics solved with forward substitution', ...
    'Forward dynamics solved with ABA'}, 'Interpreter', 'latex', 'FontSize', fontSize,  'Location','northWest','Orientation','vertical');
ylabel({'Number of floating point multiplications and additions'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'DefaultAxesFontName', 'Times New Roman');
xlabel('$N_B$ (number of links)', 'Interpreter', 'latex', 'FontSize', fontSize);

save2pdf('luAB.pdf',fH,600);

fH = figure;
fontSize = 12;
axes('Units', 'normalized', 'Parent',fH, 'FontSize', fontSize);
hold on;
plot([t_LUABA t_PLUQABA t_matrixABA ], 'LineWidth', 6)
grid on;
legend({'Forward dynamics solved with LU factorization', ...
    'Forward dynamics solved with permuted LU factorization', ...
    'Forward dynamics solved with forward substitution'}, 'Interpreter', 'latex', 'FontSize', fontSize,  'Location','northWest','Orientation','vertical');
ylabel({'Execution time [sec]'}, 'Interpreter', 'latex', 'FontSize', fontSize, 'DefaultAxesFontName', 'Times New Roman');
xlabel('$N_B$ (number of links)', 'Interpreter', 'latex', 'FontSize', fontSize);

save2pdf('t_luAB.pdf',fH,600);
