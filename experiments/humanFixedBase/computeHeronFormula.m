function [HeronFormula] = computeHeronFormula( a,b,c)
%COMPUTEHERONFORMULA  Computes Heron Formula given triangle sides having 
%lengths a, b, and c

semiper = (a + b + c)./2;
HeronFormula = sqrt(semiper*(semiper-a)*(semiper-b)*(semiper-c));

end