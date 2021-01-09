% r_VDV_compare.m% 22 July 2011%  2 March 2018 -- revised for nonlinear process% 29 March 2018 -- compare nonlinear and linear plant responses (using %                  linear MPC)% van de vuuse reactor model% first, setpoint change only  ysp = input('Setpoint change [0.1] = ') % deviation variable%  ysp = [0.17];  % setpoint change (from 0 vector)  timesp = 1;   % time of setpoint change  dist = [0];   % magnitude of input disturbance; dimension nd  timedis = 6;  % time of input disturbance  delt = 0.1;              % sample time  tfinal = input('Final time [6] = ')%  tfinal = 15;% model parameters (continuous time, state space Van de Vusse)  am = [-2.4048 0;0.8333 -2.2381];  bm = [7 7; -1.117 -1.117]; % 1 manip input, 1 dist  cm = [0 1];  dm = [0 0];  ninputs = size(bm,2);     % includes disturbance input  noutputs = size(cm,1);  nstates = size(am,1);     % number of model states  sysc_mod = ss(am,bm,cm,dm);%% discretize the model  sysd_mod = c2d(sysc_mod,delt);  [phi_mod,gamma_stuff,cd_mod,dd_stuff] = ssdata(sysd_mod)%  gamma_mod  = gamma_stuff(:,1);   % first input is manipulated  gammad_mod = gamma_stuff(:,2);   % second disturbances% ------------ plant parameters -------------------------------  planteqns = 'NL_VDVode'    % nonlinear simulation first  cplant = [0 1];   % state 2 is the output  nstates_p = 2;    % two plant states   parvec(1) = 0.5714 ; % dilution rate, 0.5714 min^-1   parvec(2) = 10;  % feed concentration of A, 10 gmol/liter   parvec(3) = 3;   % conc A, 3gmol/liter   parvec(4) = 1.117; % conc B, 1.117 gmol/liter   parvec(5) = 5/6; % k1 = 5/6 min^-1   parvec(6) = 5/3; % k2 = 5/3   parvec(7) = 1/6; % k3 = 1/6% controller parameters  p = 10;       % prediction horizon  m = 3;        % control horizon  ny = 1;       % number of measured outputs  nu = 1;       % number of manipulated inputs  nd = 1;       % number of actual disturbances  nd_est = 1;   % number of estimated disturbances  weightu = [0]; % weighting matrix for control action  weighty = [1]; % weighting matrix for outputs% constraints  umin = [-1000];  umax = [1000];  dumin = [-1000];  dumax = [1000];% Kalman Filter matrices  Q  = eye(nd_est,nd_est);       % state covariance - need to generalize!  R  = eye(ny);%    isim = 1; % additive output disturbance  iqp  = 1; % unconstrained solution  noisemag = zeros(ny,1); % no noise% run simulation  QSSmpcNLPlant% zero-order hold on input and setpoint  [tt,uu] = stairs(t,u');  [ttr,rr] = stairs(t,r');%  plot the actual plant output without measurement noise  Nfig = input('Number of first figure [1] = ')  figure(Nfig)  subplot(2,1,1)  plot(ttr,rr+1.117,':',t,y+1.117)  ylabel('y')  xlabel('time')  title(['NL VDV. DMC: P = ',num2str(p),' M = ',num2str(m)])  subplot(2,1,2)  plot(tt,uu+0.5714)  ylabel('u')  xlabel('time')  t1 = t; y1 = y+1.117; tt1 = tt; uu1 = uu+0.5714; ttr1 = ttr; rr1 = rr+1.117;% ---------------- linear simulation -----------------  planteqns = 'linVDVode'    % linear plant simulation  parvec = []; % parameters embedded in ODE file  % run simulation  QSSmpcNLPlant% zero-order hold on input and setpoint  [tt,uu] = stairs(t,u');  [ttr,rr] = stairs(t,r');  figure(Nfig+1)  subplot(2,1,1)  plot(t,y+1.117,ttr,rr+1.117,':')  ylabel('y')  xlabel('time')  title(['Linear VDV. DMC: P = ',num2str(p),' M = ',num2str(m)])  subplot(2,1,2)  plot(tt,uu+0.5714)  ylabel('u')  xlabel('time')  t2 = t; y2 = y+1.117; tt2 = tt; uu2 = uu+0.5714; ttr2 = ttr; rr2 = rr+1.117;% --------------- compare  figure(Nfig+2)  subplot(2,1,1)  plot(t1,y1,t2,y2,ttr2,rr2,':')  legend('NL','Linear')  ylabel('y')  xlabel('time')  title(['NL vs. Linear VDV plant. P = ',num2str(p),' M = ',num2str(m)])  subplot(2,1,2)  plot(tt1,uu1,tt2,uu2)  legend('NL','Linear')  ylabel('u')  xlabel('time')