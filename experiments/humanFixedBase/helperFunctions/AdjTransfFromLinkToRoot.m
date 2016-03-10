function root_X_linkNum = AdjTransfFromLinkToRoot (dmodel, q, linkNum)

% AdjTransfFromLinkToRoot provides the adjoint matrix X between the root
% reference frame and the reference frame associated to the i-th link in 
% the model.  It works both for chains and trees in Featherstone-like 
% notation.

% ASSUMPTIONS:
% - linkNum is the index (and not an ID);
% - 0 is the root;
% - Featherstone code: i_X_lamba_i = i_X_(lambda_i,i)*(lambda_i,i)_X_i = Xj*Xtree

% GOAL: 
% obtain  root_X_linkNum = 0_X_1 * 1_X_2 * ... * (lambda_linkNum)_X_linkNum


root_X_linkNum = eye(6);
i = linkNum;

if linkNum > dmodel.NB
   disp ( 'Error: choose a value of linkNum between 0 and NB!' );
   return;
end

    while (i ~= 0)
        % using jcalc for obtaining Xj 
        [Xj,~] = jcalc( dmodel.jtype{i}, q(i));

        i_X_lambda_i = Xj * dmodel.Xtree{i}; 

        % inverting i_X_lambda_i to lambda_i_X_i:
        lambda_i_X_i = InverseAdjTransform(i_X_lambda_i);

        root_X_linkNum = lambda_i_X_i * root_X_linkNum;
        i = dmodel.parent(i);
    end
end

%% Code creating for kinematic chain (not for a tree)

% % Featherstone code: 
% % i_X_lamba_i = i_X_(lambda_i,i)*(lambda_i,i)_X_i = Xj*Xtree
% 
% for i = 1 : linkNum
%     
%     % using jcalc for obtaining Xj 
%     [Xj,~] = jcalc( dmodel.jtype{i}, q(i));
%  
%     i_X_lambda_i = Xj * dmodel.Xtree{i}; 
%    
%     % inverting i_X_lambda_i to lambda_i_X_i:
%     lambda_i_X_i = zeros(6,6);
% 
%     lambda_i_X_i(1:3,1:3) = i_X_lambda_i(1:3,1:3)'; %R'
%     lambda_i_X_i(4:6,4:6) = i_X_lambda_i(4:6,4:6)'; %R'
%     lambda_i_X_i(4:6,1:3) = i_X_lambda_i(4:6,1:3)'; %S(r)R'
%  
%     root_X_linkNum = root_X_linkNum * lambda_i_X_i;
%     
% end

