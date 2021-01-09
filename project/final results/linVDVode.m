function xdot = linVDVode(t,x,flag,parvec,u,d)
%
% b.w. bequette - 22 July 2011
% linear VDV reactor model
  a = [-2.4048 0;0.8333 -2.2381];
  b = [7 7; -1.117 -1.117]; % 1 manip input, 1 dist
%
  xdot = a*x + b*[u(1);d(1)];
  