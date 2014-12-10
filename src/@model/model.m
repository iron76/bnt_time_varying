% model is a class for representing articulated rigid body models.
%
% model is a class to represent arbitrary articulated rigid body models. At
% present the class wraps the model used in Featherstone's toolbox
% (http://royfeatherstone.org/spatial/v2/download.html). The class only
% contains two properties corresponding to the number of degrees of freedom
% and the fatherstone's model
%
% PROPERTIES
%    n           - the number of degrees of freedom
%    modelParams - the Featherstone's model (as defined in the toolbox)
%
% Author: Francesco Nori
% Genova, Dec 2014


classdef model
   properties(SetAccess = immutable, GetAccess = public)
      n, S
   end
   
   properties(SetAccess = private, GetAccess = public)
      modelParams
   end
   
   methods
      function a = model(modelParams)
         % deterministicIDsolver Constructor function
         if nargin == 1
            if ~checkModel(modelParams)
               error('You should provide a featherstone-like mdoel')
            end
            a.modelParams = modelParams;
            a.n      = modelParams.NB;
            a.S      = zeros(6, a.n);

            for i = 1:a.n
               a.S(:, i) = scalc( modelParams.jtype{i} );
            end
         else
            error('You should provide a featherstone-like model')
         end
      end
   end
end
