function root_XStar_linkNum = AdjTransfStarfFromLinkToRoot (dmodel, q, linkNum)

% AdjTransfStarFromLinkToRoot provides the adjoint matrix XStar between the root
% reference frame and the reference frame associated to the i-th link in 
% the model.  It works both for chains and trees in Featherstone-like 
% notation.

% ASSUMPTIONS:
% - linkNum is the index (and not an ID);
% - 0 is the root;
% - Featherstone code: i_X_lamba_i = i_X_(lambda_i,i)*(lambda_i,i)_X_i = Xj*Xtree

% GOAL: 
% if root_X_link = [A B; C D]
% obtain  root_XStar_linkNum = [A C;B D] 

%%
root_X_linkNum = AdjTransfFromLinkToRoot(dmodel,q,linkNum);
root_XStar_linkNum = [root_X_linkNum(1:3,1:3) root_X_linkNum(4:6,1:3);...
                      root_X_linkNum(1:3,4:6) root_X_linkNum(4:6,4:6)];

end
