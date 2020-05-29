
clear all
close all

% An implementation of direct collocation
% Joris Gillis, 2018
addpath /home/thomas/Matlab/CasADI/casadi-linux-matlabR2014b-v3.5.1
addpath ./HeatFunctions
import casadi.*

% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau = collocation_points(d, 'legendre');

% Collocation linear maps
[C,D,B] = collocation_coeff(tau);
   %######## What are C,D and B ? Legendre time found ?  ######### 

% Time horizon
T = 10;

% Declare model variables
dotW = SX.sym('dotW');
Teout = SX.sym('Teout');
Tcout = SX.sym('Tcout');
Tcin = SX.sym('Tcin');
Troom = SX.sym('Troom');
x = [Teout;Tcout;Troom];
dotW = SX.sym('dotW');

% Constant parameters of the equations
Tlake= 7;
dm_w1 = 0.05; % Massflow rate of refrigerant in condenser (kg/s)
C_e = 4180 ; % Heat capacity of the evaporator (J/kg K)
m_e = 1 ; % mass of evaporator
dm_w2 = 0.05 ; % Massflow rate of refrigerant in condenser (kg/s)
C_c = 4180 ; % Heat capacity of the Condenser (J/kg K)
m_c = 1 ; % mass of condenser
C_f = 1.107 ; % Heat capacity of the refrigerant (J/kg K)
dW = 5 ; % dif of the work done to the HP cycle -> Energy input in Joules
m_water = dm_w2*T;
Cwater = 4180;
CairVbuilding = 10000;
rhoBuilding = 1.225; % kg/m³
Ke = dm_w1*C_f/C_e/m_e;
Kc = dm_w2*C_f/C_c/m_c;
Hd = 411.910431315784;
dotHd = Hd/24; % Energy divided per hours on one day !
Troom = 21;
rho_air = 1.225;
K = 1;

% Model equations

COP = (K*(0.5*(Tcin+Tcout)+273.15)/(0.5*(Tcin+Tcout)-0.5*(Tlake+Teout)));

xdot = [Ke*(Tlake-Teout)-dotW/(C_e*m_e)*((K*(0.5*(Tcin+Tcout)+273.15)/(0.5*(Tcin+Tcout)-0.5*(Tlake+Teout)))-1); ...
    Kc*(-Tcout+Tcin) + dotW/(C_c*m_c)*(K*(0.5*(Tcout+Tcin)+273.15)/(0.5*(Tcout+Tcin)-0.5*(Teout+Tlake))-1); ...
    ((Tcin-Troom)*dm_w2*Cwater-dotHd)/(CairVbuilding*rho_air)];

% Objective term
L = COP*dotW;

% Continuous time dynamics
f = Function('f', {x, dotW}, {xdot, L}) %x = [Teout;Tcout;Troom];

% Control discretization
N = 20; % number of control intervals --> How many should I choose ??  
h = T/N;

% Start with an empty NLP

opti = Opti();
J = 0;

% "Lift" initial conditions
Xk = opti.variable(3);  % Xk = [Teout, Tcout, Troom]
opti.subject_to(Xk==[10;21.1;15]);
%opti.subject_to(Tcin==26.5);

opti.set_initial(Xk, [10;21.1;15]); % not initial value, initial guess

% Collect all states/controls
Xs = {Xk};
Ws = {};


% Formulate the NLP
for k=0:N-1
   % New NLP variable for the control
   Wk = opti.variable();
   Ws{end+1} = Wk;
   opti.subject_to(0<=Wk<=50);
   %opti.subject_to(-1<=Wk<=1); % Condition de max pour Wk ? 
   %######## Ajouter contraintes pour Wk ! ######### 
      
   opti.set_initial(Wk, 0); % do we have to do this every turn ? 
   % should I put a minimum time of turning on ? 
   
   % Decision variables for helper states at each collocation point
   Xc = opti.variable(3, d);
%    opti.subject_to(-0.25 <= Xc(1,:));
   opti.set_initial(Xc, repmat([10;21.1;15],1,d));
   %######## Ajouter contraintes pour les températures !!! ######### 

   % Evaluate ODE right-hand-side at all helper states 
   [ode, Lresult] = f(Xc, Wk);

   % Add contribution to quadrature function %
   J = J + Lresult*B*h;
    % B Comes from collocation_coeff(tau) ---    %######## What is B?  ######### 
    % h =T/N = Time horizon/number of control intervals = time step for
    % integration ! 
   % Get interpolating points of collocation polynomial
   Z = [Xk Xc];

   % Get slope of interpolating polynomial (normalized)
   Pidot = Z*C;
   %######## What is C?  ######### 
   opti.subject_to(Pidot == h*ode);
  

   % State at end of collocation interval
   Xk_end = Z*D;
    %######## What is D ?  ######### 

  % New decision variable for state at end of interval
   Xk = opti.variable(3);
   opti.subject_to(Xk==[47;21.01;21]);
   opti.set_initial(Xk, [47;21.01;21]);
    %######## What is are these initial states for ??  ######### 
   Xk = opti.variable(3);
   Xs{end+1} = Xk; % add into the states variables Xs for saving! 
%    opti.subject_to(-0.25 <= Xk(1));
   opti.set_initial(Xk, [0;0;0]);
    %######## What is are these constraints for ?? for continuity ? or should I repeat every timestep the constraints ?  ######### 

   
   % Continuity constraints
   opti.subject_to(Xk_end==Xk)
   %######## What is this constraint ? What is Xk_end ?   ######### 
end

% opti.subject_to(-0.25 <= Xk(1)); T.D: For example

   %######## Add constraints !!    ######### 


Xs = [Xs{:}];
Ws = [Ws{:}];

opti.minimize(J);
opti.solver('ipopt');

sol = opti.solve();

x_opt = sol.value(Xs);

W_opt = sol.value(Ws);

% Plot the solution
% tgrid = linspace(0, T, N+1);
% clf;
% hold on
% plot(tgrid, x_opt(1,:), '--')
% plot(tgrid, x_opt(2,:), '-')
% plot(tgrid, x_opt(3,:), '-')
% plot(tgrid, x_opt(4,:), '-')
% stairs(tgrid, [W_opt nan], '-.')
% xlabel('t')
% legend('Textract';'Theatout','Theatin','Troom',W)
