function [ sigma_ygivend, logL_EM ] = computeEM(map, q, dq , Y, y, b_y, sigma_ygivend, varargin)
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

%% 
options = struct(   ...
   'MAX_ITER',      200,... maximum number of iterations
   'EPS',         1e-3);   % Convercence tolerance

%# read the acceptable names
optionNames = fieldnames(options);

%# count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('COMPUTEEM needs propertyName/propertyValue pairs')
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

%%
% Set generic parameters
len = length(q);
n_y = length(sigma_ygivend);


%% EM
% Set the initial condition for sigma_y|d
sigma_ygivend_old = sigma_ygivend;

Sigma_EM = cell(options.MAX_ITER,1);  
% Loop for EM iterations
for n = 1 : options.MAX_ITER

    A = zeros(size(sigma_ygivend_old));  
    for i = 1:len    
        
        % compute MAP
        [mu_dgiveny, Sigma_dgiveny, map] = computeMAP(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), sigma_ygivend_old);
        yHat = y(:,i) - b_y(:,i);
        A = A + (Y{i}*Sigma_dgiveny*Y{i}' + (yHat-(Y{i}*mu_dgiveny))*(yHat-(Y{i}*mu_dgiveny))');

    end
    sigma_ygivend = (1/len)*A;   
        
    

      % exit condition
       fprintf('EM iteration %i :norm %i\n',n, norm(sigma_ygivend - sigma_ygivend_old));
      if abs(norm(sigma_ygivend - sigma_ygivend_old)) < options.EPS;
         fprintf('Local convergence reached with %i iterations.', n);
         break
      else
          sigma_ygivend_old = sigma_ygivend;
          Sigma_EM{n} = sigma_ygivend;
      end        
end

%% LogLikelihood

logL_EM = zeros(options.MAX_ITER,1);
for k =1:options.MAX_ITER

% % Set the initial condition for the log-likelihood
% logL_old = -inf;
% % for i = 1:len  
% %     logL_old = logL_old +  0.5*( n_y*log(2*pi) + log(det(inv(sigma_ygivend_old))) - (y(:,i))'*sigma_ygivend_old*(y(:,i)) );
% % end

    logL = 0;
    for i = 1:len 
        
        [~, ~, map] = computeMAP(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), Sigma_EM{k});
        yHat = y(:,i) - b_y(:,i);
        % compute LogL
        logL = logL +  ( -n_y*log(2*pi) + log(det(inv(Sigma_EM{k} + Y{i}*map.sigmaBarD*Y{i}'))) - (yHat-Y{i}*map.muBarD)'*inv(Sigma_EM{k} + Y{i}*map.sigmaBarD*Y{i}')* (yHat-Y{i}*map.muBarD));
    end
    logL_EM(k,1)= 0.5*logL;
    
end

%%
% % % % Set generic parameters
% % % len = length(q);
% % % sizeD = 0;
% % % sizeY = 0;
% % % 
% % % % 1. set the initial condition for sigma_y|d
% % % sigma_ygivend_old = sigma_ygivend;
% % % logL_old = 0;
% % % 
% % % for i = 1:len   
% % %     [mu_dgiveny, ~, map] = computeMAP(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), sigma_ygivend_old);
% % %     yHat = y(:,i) - b_y(:,i);
% % %     logL_old = logL_old + ( -(52+26)*log(2*pi) + log(det(full(inv(map.sigmaBarD)))) + log(det(inv(sigma_ygivend_old))) - (mu_dgiveny-map.muBarD)'*full(inv(map.sigmaBarD))*(mu_dgiveny-map.muBarD) - (yHat-(Y{i}*mu_dgiveny))'*inv(sigma_ygivend_old)*(yHat-(Y{i}*mu_dgiveny)));
% % %     
% % % end
% % % logL_old = 0.5*logL_old;
% % % 
% % % 
% % % fprintf('Initial log-likelihood => %f\n', logL_old);
% % % 
% % % % 2. loop for EM iterations
% % % for n = 1 : options.MAX_ITER
% % % 
% % %     A = zeros(size(sigma_ygivend_old));
% % %     
% % %     logL = 0;
% % %     
% % %     for i = 1:len    
% % %         % compute MAP
% % %         [mu_dgiveny, Sigma_dgiveny, ~] = computeMAP(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), sigma_ygivend_old);
% % %         yHat = y(:,i) - b_y(:,i);
% % %         A = A + (Y{i}*Sigma_dgiveny*Y{i}' + (yHat-(Y{i}*mu_dgiveny))*(yHat-(Y{i}*mu_dgiveny))');
% % %     end
% % %     
% % %       sigma_ygivend = (1/len)*A;
% % %       % compute log-likelihood
% % %       for i = 1:len
% % %         [mu_dgiveny, ~, map] = computeMAP(map,q(:,i),dq(:,i),Y{i}, y(:,i),b_y(:,i), sigma_ygivend);
% % %         yHat = y(:,i) - b_y(:,i);
% % %          logL = logL + 0.5*( -(52+26)*log(2*pi) + log(det(full(inv(map.sigmaBarD)))) + log(det(inv(sigma_ygivend))) - (mu_dgiveny-map.muBarD)'*full(inv(map.sigmaBarD))*(mu_dgiveny-map.muBarD) - (yHat-(Y{i}*mu_dgiveny))'*inv(sigma_ygivend)*(yHat-(Y{i}*mu_dgiveny)));
% % %       end 
% % %       
% % %       
% % % 
% % %       % 3. CHECK CONVERGENCE 
% % %       fprintf('EM iteration %i: log-likelihood => %f\n', n, logL);
% % %       if abs(logL - logL_old) < options.EPS;
% % %          fprintf('Local convergence reached with %i iterations.', n);
% % %          break
% % %       else
% % %           sigma_ygivend_old = sigma_ygivend;
% % %           logL_old = logL;
% % %       end  
% % %         
% % % end

end

