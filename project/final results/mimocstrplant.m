% Function file for Nonlinear cstr plant

function dedt = mimocstrplant(t,x,u)

% parameter vectors
 % reactor volume m3
Tf = 350; % feed temperature K
lambda = -69.71*10^6; % heat of reaction J/Kmol
Cao = 8.01; % Feed conc Kmol 
ko = 20.75*10^6; % pre exponential factor s-1
E = 69.71*10^6; % activation energy L/kmol
R = 8314;
rho = 801;
Cp = 3137;
U = 851;
Aj = 101;
Vj = 10.1;
Tcin = 294;
rho_j = 1000;
Cj = 4183;
h_ss = 1.0 ; % reactor steady state height
Ca_ss = 0.0842;
Tr_ss = 339.7022;
Tj_ss = 323.7669;
F_ss = 0.04377;
Fj_ss = 0.011;
alpha = 0.04377;

h_dev = x(1);
Ca_dev = x(2);
Tr_dev = x(3);
Tj_dev = x(4);

% steady state plus deviation for all variables 
h = h_ss + h_dev;
Ca = Ca_ss + Ca_dev;
Tr = Tr_dev + Tr_ss;
Tj = Tj_dev + Tj_ss;

F = u(1) + F_ss;
Fj = u(2) + Fj_ss;


% Ca cannot be less than zeros
 if Ca <=0
       Ca = 0;
 end
 
 % Plant equations
dhdt = (F - alpha*h)/Aj; 
dCdt = (F*Cao - alpha*h*Ca)/(Ar*h) -(Ca*(F - alpha*h))/(Aj*h) - Aj*h*Ca*ko*exp((-E)/(R*Tr));
dTrdt = (F*Tf- alpha*h*Tr)/(Aj*h) - (Tr*(F - alpha*h))/(h*Aj) - (lambda*Ca*ko*exp((-E)/(R*Tr)))/(rho*Cp) - (U*Aj*(Tr-Tj))/(Aj*h*rho*Cp);
dTjdt = (Fj/Vj)*(Tcin-Tj) + (U*Aj*(Tr-Tj))/(Vj*rho_j*Cj);
dedt = [dhdt;dCdt;dTrdt;dTjdt];
end