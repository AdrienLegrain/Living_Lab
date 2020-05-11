%
%     This file is part of CasADi.
%
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
%                             K.U. Leuven. All rights reserved.
%     Copyright (C) 2011-2014 Greg Horn
%
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
%
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
%
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%

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
x = [Textract;Theatout];

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
Ke = dm_w1*C_f/C_e/m_e;
Kc = dm_w2*C_f/C_c/m_c;
Hd = 411.910431315784;
Troom = 21;

% Model equations

xdot = [-Ke*Textract-Ke*Tlake - dotW/(C_e*m_e)*(Theatout+Theatin+2*273.15)/((Theatout+Theatin-Tlake-Textract)-1);
    -Kc*(-Theatout+Theatin) - dotW/(C_c*m_c)*(K*(0.5*(Theatout+Theatin)+273.15)/(0.5*(Theatout+Theatin)-0.5*(Textract+Tlake))-1)]
% dotTextract = -Ke*Textract-Ke*Tlake - dotW/(C_e*m_e)*(Theatout+Theatin+2*273.15)/((Theatout+Theatin-Tlake-Textract)-1);
% dotTheatout = -Kc*(-Theatout+Theatin) - dotW/(C_c*m_c)*(K*(0.5*(Theatout+Theatin)+273.15)/(0.5*(Theatout+Theatin)-0.5*(Textract+Tlake))-1);
% Troom = (Theatout-Theatin)*m_water*Cwater-Hd)/(CairVbuilding+m_water*Cwater);
% Theatin = (Troom*CairVbuilding + Hd+Cwater*m_water*Theatout)/(Cwater*m_water)
% COP = K*(0.5*(Theatout+Theatin)+273.15)/(0.5*(Theatout+Theatin)-0.5*(Textract+Tlake))

% Objective term
L = COP*dotW;

% Continuous time dynamics
f = Function('f', {x, dotW, Theatin}, {xdot, L, (Troom*CairVbuilding + Hd+Cwater*m_water*Theatout)/(Cwater*m_water)});

% --------------------VVV FUNCTION STILL NOT CHANGED HEREUNDRER VVV-------------------Thomas

% Control discretization
N = 20; % number of control intervals --> How many should I choose ??  
h = T/N;

% Start with an empty NLP

opti = Opti();
J = 0;

% "Lift" initial conditions
Xk = opti.variable(2);
opti.subject_to(Xk==[0; 1]); % what is this ? The initial condition as well ? 
opti.set_initial(Xk, [0; 1]);

% Collect all states/controls
Xs = {Xk};
Us = {};

% Formulate the NLP
for k=0:N-1
   % New NLP variable for the control
   Uk = opti.variable();
   Us{end+1} = Uk;
   opti.subject_to(-1<=Uk<=1);
   opti.set_initial(Uk, 0);

   % Decision variables for helper states at each collocation point
   Xc = opti.variable(2, d);
   opti.subject_to(-0.25 <= Xc(1,:));
   opti.set_initial(Xc, repmat([0;0],1,d));

   % Evaluate ODE right-hand-side at all helper states --> This is solving
   % the ODE ??
   [ode, quad] = f(Xc, Uk);

   % Add contribution to quadrature function % I should add the total
   % result right ? so COP*dotW ? 
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
   Xk = opti.variable(2);
   Xs{end+1} = Xk;
   opti.subject_to(-0.25 <= Xk(1));
   opti.set_initial(Xk, [0;0]);

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
tgrid = linspace(0, T, N+1);
clf;
hold on
plot(tgrid, x_opt(1,:), '--')
plot(tgrid, x_opt(2,:), '-')
stairs(tgrid, [u_opt nan], '-.')
xlabel('t')
legend('x1','x2','u')
