% MAP solves inverse dynamics with redundant measurements.
%
% The class is instantiated with two additional classes to describe the
% articulated rigid body model, e.g.:
%            myModel = model(autoTree(NB)) 
% and to describe the sensor distribution, e.g.: 
%            mySens  = sensors(autoSensStochastic(autoSensMAP(myModel)))
% with the following instanatiation: 
%             myMAP = MAP(myModel, mySens).
% Computations are then performed with the method d = myMAP.solveID();
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef dMAP < derivativeMAPsolver
   properties
      sparsified = 0;
      S
   end
   
   methods
      function b = dMAP(m,y)
         if nargin == 0
            error(['You should provide a ' ...
               'model to instantiate MAP'] )
         else
            if ~checkModel(m.modelParams)
               error(['You should provide a featherstone-like ' ...
                  'model to instantiate MAP'] )
            end
         end
         b = b@derivativeMAPsolver(m,y);                 

      end % MAP
            
      function disp(b)
         % Display MAP object
         disp@derivativeMAPsolver(b)
         fprintf('MAP disp to be implemented! \n')
      end % disp
     
   end % methods
end % classdef

