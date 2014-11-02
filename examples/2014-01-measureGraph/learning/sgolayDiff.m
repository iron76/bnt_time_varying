function [SG0, SG1 ] = sgolayDiff(y, x)
%SGOLAYDIFF differentiate with Savitzky-Golay
%   Detailed explanation goes here

N = 3;                 % Order of polynomial fit
F = 501;               % Window length
[b,g] = sgolay(N,F);   % Calculate S-G coefficients

dx = x(2)-x(1);

HalfWin  = ((F+1)/2) -1;
for n = (F+1)/2:length(x)-(F+1)/2,
  % Zero-th derivative (smoothing only)
  SG0(n) =   dot(g(:,1), y(n - HalfWin: n + HalfWin));
  
  % 1st differential
  SG1(n) =   dot(g(:,2), y(n - HalfWin: n + HalfWin));
end

SG0(length(x)-(F+1)/2+1:length(x))=0;
SG1(length(x)-(F+1)/2+1:length(x))=0;

SG1 = SG1/dx;         % Turn differential into derivative

end

