% DNEA solves inverse dynamics with redundant measurements.
%
% The class is instantiated with two additional classes to describe the
% articulated rigid body model, e.g.:
%            myModel = model(autoTree(NB)) 
% and to describe the sensor distribution, e.g.: 
%            mySens  = sensors(autoSensStochastic(autoSensSNEA(myModel)))
% with the following instanatiation: 
%             myDNEA = DNEA(myModel, mySens).
% Computations are then performed with the method d = myDNEA.solveID();
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef DNEA < derivativeIDsolver
      
   methods
      function b = DNEA(m,y)
         if nargin == 0
            error(['You should provide a ' ...
               'model to instantiate SNEA'] )
         else
            if ~checkModel(m.modelParams)
               error(['You should provide a featherstone-like ' ...
                  'model to instantiate SNEA'] )
            end
         end
         b    = b@derivativeIDsolver(m,y);
         b.Sd = zeros(b.IDmodel.modelParams.NB*26, b.IDmodel.modelParams.NB*26);
         
      end % DNEA
      
      function disp(b)
         % Display DNEA object
         disp@deterministicIDsolver(b)
         fprintf('DNEA disp to be implemented! \n')
      end % disp
      
   end % methods
end % classdef

