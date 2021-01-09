function xdot = NL_VDVode(t,x,flag,parvec,u,d)
%
% b.w. bequette - 2 March 2018
% nonlinear VDV reactor model
%  variables passed in as deviation variables
%  parvec contains the steady-state values in addition to the 
%    kinetic parameters
   Ca_dev= x(1);
   Cb_dev= x(2);
%
   FOVss = parvec(1); % dilution rate, 0.5714 min^-1
   Cafss = parvec(2); % feed concentration of A, 10 gmol/liter
   Cass  = parvec(3); % conc A, 3gmol/liter
   Cbss  = parvec(4); % conc B, 1.117 gmol/liter
   k1    = parvec(5); % 5/6 min^-1
   k2    = parvec(6); % 5/3
   k3    = parvec(7); % 1/6
%
   FOV   = FOVss + u(1); % dilution rate
   Caf   = Cafss + d(1); % feed concentration
%
   Ca    = Cass + Ca_dev;
   Cb    = Cbss + Cb_dev;
   if Ca <=0;
       Ca = 0;
   end
       if Cb <= 0;
           Cb = 0
       end
%
   dCadt = FOV*(Caf - Ca) - k1*Ca -k3*Ca*Ca;
   dCbdt = -FOV*Cb + k1*Ca - k2*Cb;%
    xdot = [dCadt;dCbdt];
  