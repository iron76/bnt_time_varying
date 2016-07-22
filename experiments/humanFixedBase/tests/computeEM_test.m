% function [Sigma_EM, l_EM] = computeEM_test(map, q, dq , Y, y, b_y, sigma_ygivend, varargin)
%function [Sigma_EM, l_EM] = computeEM_test(map, q, dq , Y, y, b_y, sigma_ygivend, varargin)
%function [ sigma_ygivend, Sigma_EM, logL_EM, Q_EM] = computeEM_test(map, q, dq , Y, y, b_y, sigma_ygivend,  sigma_Dmodel,varargin)
function [Sigma_EM, Qin_EM, lin_EM, Qfin_EM, lfin_EM] = computeEM_test(map, q, dq , Y, y, b_y, sigma_ygivend, varargin)
%% EM
% Function for computing  Expectation-Maximization algorithm. As in this 
% code it is used for the data-driven covariance estimation, thus the 
% parameters to be estimated is the covariance of measurements given vector
% d (Sigma_y|d).

% Inputs: 
% -  map                map object;
% -  q                  joint position;
% -  dq                 joint velocity;
% -  maxIterNum         maximum nuber of EM iterations;
% -  len                number of samples;
% -  epsilon            quantity for convergence condition;
% ...
% - sigma_ygivend       initial value for sigma_ygivend, i.e. the sigma
%                       coming from datasheet.

%% initialOptions
options = struct(   ...
   'MAX_ITER',      200,... maximum number of iterations
   'EPS',         1e-4);   % Convercence tolerance

%# read the acceptable names
optionNames = fieldnames(options);

%# count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('COMPUTE EM needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = upper(pair{1}); %# make case insensitive

if any(strcmp(inpName,optionNames))
   %# overwrite options. If you want you can test for the right class here
   %# Also, if you find out that there is an option you keep getting wrong,
   %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
   options.(inpName) = pair{2};
else
   error('%s is not a recognized parameter name',inpName)
   end
end

%% testOptions
ThetaSigma_ygivend = true;
ThetaSigma_Dmodel = false;

%% ========================ESTIMATION wrt theta1========================== %%

 %% ========================ESTIMATION wrt theta1========================== %%
if(ThetaSigma_ygivend)
    options.MAX_ITER = 10;
     %%EM (by using analitical derivation)

    % Set generic parameters
    len = length(q);
%     n_y = map.IDsens.m;
%     n_D = length(map.d);

    % Set the initial condition for sigma_y|d
    
    sigma_ygivend_old = sigma_ygivend;

    Sigma_EM = cell(options.MAX_ITER,1); 
    Qin_EM = zeros(options.MAX_ITER,1);
    lin_EM = zeros(options.MAX_ITER,1);
    Qfin_EM = zeros(options.MAX_ITER,1);
    lfin_EM = zeros(options.MAX_ITER,1);
    
    %     Q_paramTest = cell(options.MAX_ITER,1);
    
    % Loop for EM iterations
    for k = 1 : options.MAX_ITER

        A = zeros(size(sigma_ygivend_old));
        len = 1;
        
        % CHECK INEQUALITY PRE-EM
        Qin = 0;
        for i = 1:len    
              Q_fun = computeQ(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), sigma_ygivend_old,sigma_ygivend_old);
              Qin = Qin + Q_fun;
        end
        Qin_EM(k,1) = Qin;
        
        lin1 = 0;
        %lin2 = 0;
        for i = 1:len    
           [l_fun1] = computeLogL(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), sigma_ygivend_old,sigma_ygivend_old);
           lin1 = lin1 + l_fun1;
%            lin2 = lin2 + l_fun2;
        end 
        lin_EM(k,1) = lin1;
%         lin_EM(k,2) = lin2;

        % =============compute Sigma
        for i = 1:len    
            % compute MAP
            [mu_dgiveny, Sigma_dgiveny, map] = computeMAP(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), inv(sigma_ygivend_old));
            yHat = y(:,i) - b_y(:,i);
            
            A1 = 2 * ((Y{i}*mu_dgiveny) - yHat)*((Y{i}*mu_dgiveny) - yHat)';
            A2 = -0.5 * diag(diag(A1));
            A3 = 2*Y{i}*Sigma_dgiveny*Y{i}';
            A4 = - 0.5 * diag(diag(A3));
            A = A + (A1 + A2 + A3 + A4);
        end
        sigma_ygivendBar = (1/len)*A;       
        sigma_ygivend = 0.5*sigma_ygivendBar;
        sigma_ygivend = sigma_ygivend + diag(diag(sigma_ygivend));
        Sigma_EM{k} = sigma_ygivend;
        % ==========================
        
% % %         % =============compute Sigma
% % %         for i = 1:len    
% % %             compute MAP
% % %             [mu_dgiveny, Sigma_dgiveny, map] = computeMAP(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), inv(sigma_ygivend_old));
% % %             yHat = y(:,i) - b_y(:,i);
% % %             
% % %             clatest
% % %             A1 = yHat*(yHat') - 0.5 *diag(diag(yHat*yHat'));
% % %             A2 = Y{i}*(Sigma_dgiveny+(mu_dgiveny*mu_dgiveny'))*Y{i}'- 0.5 * diag(diag(Y{i}*(Sigma_dgiveny+(mu_dgiveny*mu_dgiveny'))*Y{i}'));
% % %             A3  = -(Y{i}*mu_dgiveny*yHat') - (yHat*mu_dgiveny'*Y{i}') + (diag(diag(Y{i}*mu_dgiveny*yHat')));
% % %             A = A + (A1 + A2 + A3);
% % %         end
% % %         sigma_ygivendBar = (1/len)*A;       
% % %         sigma_ygivend = sigma_ygivendBar + diag(diag(sigma_ygivendBar));
% % %         sigma_ygivend = sigma_ygivendBar + diag(diag(sigma_ygivendBar));
% % %         sigma_ygivend = 0.5 * (sigma_ygivend + sigma_ygivend');
% % %         Sigma_EM{k} = sigma_ygivend;
% % %         ==========================
        
        % CHECK INEQUALITY POST-EM
        Qfin = 0;
        for i = 1:len    
              Q_fun = computeQ(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), sigma_ygivend_old,sigma_ygivend);
              Qfin = Qfin + Q_fun;
        end
        Qfin_EM(k,1) = Qfin;
               
        
        lfin1 = 0;
%         lfin2 = 0;
        for i = 1:len    
           [l_fun1] = computeLogL(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), sigma_ygivend,sigma_ygivend);
           lfin1 = lfin1 + l_fun1;
%            lfin2 = lfin2 + l_fun2;
        end 
        lfin_EM(k,1) = lfin1;
%         lfin_EM(k,2) = lfin2;
        
        % exit condition
        fprintf('EM iteration %i :norm %i\n',k, norm(sigma_ygivend - sigma_ygivend_old));
        if abs(norm(sigma_ygivend - sigma_ygivend_old)) < options.EPS;
           fprintf('Local convergence reached with %i iterations.', k);
           break
        else
           sigma_ygivend_old = sigma_ygivend;
        end

    end
   
%% ========================ESTIMATION wrt theta2========================== %%
if(ThetaSigma_Dmodel)
    
    %% EM (by using analitical derivation)

    % Set generic parameters
    len = length(q);
    n_y = length(sigma_ygivend);

    % Set the initial condition for sigma_Dmodel
    sigma_Dmodel_old = sigma_Dmodel;
   

    Sigma_EM = cell(options.MAX_ITER,1);  
    % Loop for EM iterations
    for n = 1 : options.MAX_ITER
        A = zeros(52);
    %sigma_ygivendBar = zeros(size(sigma_ygivend));
    %A = zeros(size(sigma_ygivend_old));
        for i = 1:len    

            % compute MAP
            [mu_dgiveny, Sigma_dgiveny, map] = computeMAP(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), inv(sigma_ygivend), inv(sigma_Dmodel_old));

            A1 = 2*(Sigma_dgiveny+(mu_dgiveny*mu_dgiveny'))-diag(diag(Sigma_dgiveny+(mu_dgiveny*mu_dgiveny')));
            A2 = - 2*map.muBarD*mu_dgiveny' - 2*mu_dgiveny*map.muBarD' + 2*diag(diag(map.muBarD*mu_dgiveny'));
            A3  = + 2*map.muBarD*map.muBarD' - diag(diag(map.muBarD*map.muBarD'));
            A = A + (A1 + A2 + A3);
        end

        sigma_Dmodel = (1/len)*A;       
        sigma_Dmodel = 0.5*sigma_Dmodel;
        sigma_Dmodel = sigma_Dmodel + diag(diag(sigma_Dmodel));
        Sigma_EM{n} = sigma_Dmodel;
        % exit condition
        fprintf('EM iteration %i :norm %i\n',n, norm(sigma_Dmodel - sigma_Dmodel_old));
        if abs(norm(sigma_Dmodel - sigma_Dmodel_old)) < options.EPS;
           fprintf('Local convergence reached with %i iterations.', n);
           break
        else
           sigma_Dmodel_old = sigma_Dmodel;
        end
        
        Sigma_EM{n} = sigma_Dmodel;
    end


end
end
