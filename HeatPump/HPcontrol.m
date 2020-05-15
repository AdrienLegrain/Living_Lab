

% An implementation of direct collocation
% Joris Gillis, 2018
import casadi.*

% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau = collocation_points(d, 'legendre');

% Collocation linear maps
[C,D,B] = collocation_coeff(tau);

% Time horizon
T = 10;

% Declare model variables
dotW = SX.sym('dotW');
Textract = SX.sym('Textract');
Theatout = SX.sym('Theatout');
Theatin = SX.sym('Theatin');
x = [Textract;Theatout,Theatin,Troom];

% Constant parameters of the equations
Tlake= 7;
dm_w1 = 11.38; % Massflow rate of refrigerant in condenser (kg/s)
C_e = 4180 ; % Heat capacity of the evaporator (J/kg K)
m_e = 1 ; % mass of evaporator
dm_w2 = 11.38 ; % Massflow rate of refrigerant in condenser (kg/s)
C_c = 4180 ; % Heat capacity of the Condenser (J/kg K)
m_c = 1 ; % mass of condenser
C_f = 1.107 ; % Heat capacity of the refrigerant (J/kg K)
dW = 5 ; % dif of the work done to the HP cycle -> Energy input in Joules
m_water = dm_w2*T;
Cwater = 4180;
CairVbuilding = 10000;
rhoBuilding = 1.225 % kg/m³
Ke = dm_w1*C_f/C_e/m_e;
Kc = dm_w2*C_f/C_c/m_c;
Hd = 411.910431315784;
dotHd = Hd/24; % Energy divided per hours on one day !
Troom = 21;
rho_air = 1.225;
K = 1;
% Model equations
xdot = [dotTextract;dotTheatout;Theatin;Troom] 
xdot = [-Ke*Textract+Ke*Tlake-dotW/(C_e*m_e)*(K*(0.5*(Theatout+Theatin)+273.15)/(0.5*(Theatout+Theatin)-0.5*(Textract+Tlake))-1);
    Kc*(-Theatout+Theatin) - dotW/(C_c*m_c)*(K*(0.5*(Theatout+Theatin)+273.15)/(0.5*(Theatout+Theatin)-0.5*(Textract+Tlake))-1);
    (Theatin-Theatout)*dm_w2*Cwater-Hd)/(CairVbuilding*rhoBuilding);
    Theatout-0.1] % CONSTANT ASSUMPTION !!! 

% dotTextract = -Ke*Textract+Ke*Tlake-dotW/(C_e*m_e)*(K*(0.5*(Theatout+Theatin)+273.15)/(0.5*(Theatout+Theatin)-0.5*(Textract+Tlake))-1)
% dotTheatout = Kc*(-Theatout+Theatin) - dotW/(C_c*m_c)*(K*(0.5*(Theatout+Theatin)+273.15)/(0.5*(Theatout+Theatin)-0.5*(Textract+Tlake))-1);
% Theatin = (CairVbuilding*rho_air*dotTroom + dotHd)/(Cwater*dmw_2)+Theatout;
% dotTroom = ((Theatin-Theatout)*Cwater*dmw_2-dotHd)/CairVbuilding*rho_air;
% COP = K*(0.5*(Theatout+Theatin)+273.15)/(0.5*(Theatout+Theatin)-0.5*(Textract+Tlake))


% Objective term
L = COP*dotW;

% Continuous time dynamics
f = Function('f', {x, dotW}, {xdot, L})

% Control discretization
N = 20; % number of control intervals --> How many should I choose ??  
h = T/N;

% Start with an empty NLP

opti = Opti();
J = 0;

% "Lift" initial conditions
Xk = opti.variable(4);
opti.subject_to(Xk==[47;21.01;52;21]);
opti.set_initial(Xk, [47;21.01;52;21]);

% Collect all states/controls
Xs = {Xk};
Ws = {};


% Formulate the NLP
for k=0:N-1
   % New NLP variable for the control
   Wk = opti.variable();
   Us{end+1} = Wk;
   opti.subject_to(-1<=Wk<=1); % Condition sur W ? 
   opti.set_initial(Wk, 0);

   % Decision variables for helper states at each collocation point
   Xc = opti.variable(4, d);
   opti.subject_to(-0.25 <= Xc(1,:));
   opti.set_initial(Xc, repmat([0;0;0;0],1,d));

   % Evaluate ODE right-hand-side at all helper states --> This is solving
   % the ODE ??
   [ode, quad] = f(Xc, Wk);

   % Add contribution to quadrature function %
   J = J + quad*B*h;

   % Get interpolating points of collocation polynomial
   Z = [Xk Xc];

   % Get slope of interpolating polynomial (normalized)
   Pidot = Z*C;
   % Match with ODE right-hand-side 
   
   opti.subject_to(Pidot == h*ode);
  

   % State at end of collocation interval
   Xk_end = Z*D;

   % New decision variable for state at end of interval
   Xk = opti.variable(4);
   opti.subject_to(Xk==[47;21.01;52;21]);
   opti.set_initial(Xk, [47;21.01;52;21]);
   
   Xk = opti.variable(4);
   Xs{end+1} = Xk;
   opti.subject_to(-0.25 <= Xk(1));
   opti.set_initial(Xk, [0;0;0;0]);

   % Continuity constraints
   opti.subject_to(Xk_end==Xk)
end


Xs = [Xs{:}];
Us = [Us{:}];

opti.minimize(J);
opti.solver('ipopt');

sol = opti.solve();

x_opt = sol.value(Xs);

u_opt = sol.value(Us);

% Plot the solution
% tgrid = linspace(0, T, N+1);
% clf;
% hold on
% plot(tgrid, x_opt(1,:), '--')
% plot(tgrid, x_opt(2,:), '-')
% stairs(tgrid, [u_opt nan], '-.')
% xlabel('t')
% legend('x1','x2','u')
