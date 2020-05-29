clear all
close all

% An implementation of direct collocation
% Joris Gillis, 2018
addpath /home/thomas/Matlab/CasADI/casadi-linux-matlabR2014b-v3.5.1
addpath ./HeatFunctions
import casadi.*





%% Ex 3 MPC problem setup

opti = casadi.Opti(); % Optimization problem

N = 10; % MPC horizon

HPdynamics
h = 0.025;
f_HP = @(x,u) RK4(x,u,h,f);
% The function f takes (X,U) as argument 

%% ---- decision variables ---------
X = opti.variable(3); % state trajectory variables
dotW = opti.variable(1);
U = opti.variable(2);   % control trajectory (Tcin, dotHd)

X0   = opti.parameter(3,1); % parameter variable for initial state
opti.set_value(X0, [10;21.1;15])
Teout0   = X0(1);
Tcout0 = X0(2); 
Troom0 = X0(3);


dotW0 = 0;

Teout   = X0(1);
Tcout = X0(2); 
Troom = X0(3);
dotW = dotW0;

% epsilon_speed = opti.variable(1,N+1); % slack variable for speed constraint
COP = K*(0.5*(U(1)+X(2))+273.15)/(0.5*(U(1)+X(2))-0.5*(Tlake+X(1)));

% ---- objective ---------
opti.minimize(COP*dotW);

% ---- multiple shooting --------
for k=1:N % loop over control intervals
  
    opti.subject_to(X(:,k+1) == f_HP(X(:,k), U(:,k)));
  
  %%%% WRITE YOUR DYNAMICS CONSTRAINT HERE
  %   opti.subject_to( ... );
  
end

% ---- path constraints -----------

% limit = track.maxspeed;
% opti.subject_to(speed  <=   limit(pos) + epsilon_speed); % track speed limit
% opti.subject_to(0 <= U <= 1);  % control is limited
% 
% % ---- boundary conditions --------
% opti.subject_to(pos(1)==pos0);   % use initial position
% opti.subject_to(speed(1)==v0); % use initial speed
% opti.subject_to(epsilon_speed >= 0); % slack lower bound
% 
% % Pass parameter values
% opti.set_value(pos0, 0.0);
% opti.set_value(v0, 0);

% ---- Setup solver NLP    ------
ops = struct;
ops.ipopt.print_level = 0;
ops.ipopt.tol = 1e-3;
opti.solver('ipopt', ops);
sol = opti.solve();   % actual solve

% % ---- Plot predicted trajectory for debugging ------
% 
% figure(1); clf
% title('Predicted trajectory from initial state')
% hold on; grid on
% plot(limit(sol.value(pos)),'k--','linewidth',2);
% plot(sol.value(pos),'linewidth',2);
% plot(sol.value(speed),'linewidth',2);
% stairs(1:N, sol.value(U(1,:)),'linewidth',2);
% stairs(1:N,-sol.value(U(2,:)),'linewidth',2);
% legend('speed limit','pos','speed','throttle','brake','Location','northwest')

%% Functions

function [x_next] = RK4(X,U,h,f)
%
% Inputs : 
%    X, U current state and input
%    h    sample period
%    f    continuous time dynamics f(x,u)
% Returns
%    State h seconds in the future
%

% Runge-Kutta 4 integration
% write your function here
   k1 = f(X,         U);
   k2 = f(X+h/2*k1, U);
   k3 = f(X+h/2*k2, U);
   k4 = f(X+h*k3,   U);
   x_next = X + h/6*(k1+2*k2+2*k3+k4);
end

% function plotSims(t, sims)
% 
% figure(1); clf
% colors = get(gca,'ColorOrder');
% markers = {'o','d','s','^','<','p','h'};
% for i = 1:length(sims)
%   sims{i}.style = {'linewidth',2,'marker',markers{i},'markersize',5,'markerfacecolor','w','color',colors(i,:)};
% end
% 
% leg = {};
% 
% for i = 1:length(sims)
%   sim = sims{i};
%   
%   try
%     subplot(2,1,1)
%     plot(t,sim.X(1,:),sim.style{:});
%     hold on; grid on
%     ylabel('\lambda')
%     
%     subplot(2,1,2)
%     plot(t,sim.X(2,:),sim.style{:});
%     hold on; grid on
%     xlabel('t (s)')
%     ylabel('v')
%     
%     leg{end+1} = sim.name;
%   catch
%     warning('Define the structure %s above\n', sim.name);
%   end
% end
% legend(leg{:})
% end
