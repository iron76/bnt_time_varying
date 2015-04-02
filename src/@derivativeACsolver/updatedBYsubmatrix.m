function [ obj ] = updatedBYsubmatrix( obj )
%UPDATESTATEDERIVATIVESUBMATRIX Compute \frac{\partial (Dd + b)}{\partial x}
%   Compute the derivative of D(q)d + b(q,\dot{q}) with respect to x
%   (x is defined as (q,\dot{q}).
%   The output is saved in the dDb (please suggest a better name) attribute
%   of the input obj.
%

for j = 1 : length(obj.idby_s)
   sens_str = obj.IDsens.sensorsParams.labels{obj.idby_s(j),1};
   h  = str2double(sens_str(8:end));
   dv = obj.dvdx{h, obj.jdby_s(j)};
   obj.dby_s = set(obj.dby_s, dv(1:3,1), obj.idby_s(j), obj.jdby_s(j));
end