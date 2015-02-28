function [ obj ] = updateBYsubmatrix( obj )
%UPDATESTATEDERIVATIVESUBMATRIX Compute \frac{\partial (Dd + b)}{\partial x}
%   Compute the derivative of D(q)d + b(q,\dot{q}) with respect to x
%   (x is defined as (q,\dot{q}).
%   The output is saved in the dDb (please suggest a better name) attribute
%   of the input obj.
%

for i = 1 : length(obj.iby_s)
   sens_str = obj.IDsens.sensorsParams.labels{obj.iby_s(i),1};
   j = str2double(sens_str(8:end));
 
   obj.by_s = set(obj.by_s, obj.v(1:3,j), obj.iby_s(i), obj.jby_s(i));
end