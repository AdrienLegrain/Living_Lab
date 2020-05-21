import casadi.*


h = 0.01; % integration time step
N = 10; % MPC horizon



%% FUNCTION DEFINITIONS
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

h = 0.01; % integration time step

f_discrete = @(x,u) RK4(x,u,h,f);

%% Ex 3 MPC problem setup

opti = casadi.Opti(); % Optimization problem

% ---- decision variables ---------
X = opti.variable(2,N+1); % state trajectory variables
U = opti.variable(2,N);   % control trajectory (throttle, brake)

X0   = opti.parameter(2,1); % parameter variable for initial state
pos0 = X0(1); % Initial position and velocity
v0   = X0(2);

pos   = X(1,:);
speed = X(2,:);

epsilon_speed = opti.variable(1,N+1); % slack variable for speed constraint

% ---- objective ---------
opti.minimize(...
  -10*sum(X(2,:))  + ... % Max velocity
  0.1*U(1,:)*U(1,:)' + ... % Minimize accel
  10*U(2,:)*U(2,:)'   + ... % Minimize braking
  10000*(epsilon_speed(1,:)*epsilon_speed(1,:)' + sum(epsilon_speed))); % Soft constraints

% ---- multiple shooting --------
for k=1:N % loop over control intervals
  
  opti.subject_to(X(:,k+1) == f_discrete(X(:,k), U(:,k)));
  
  %%%% WRITE YOUR DYNAMICS CONSTRAINT HERE
  %   opti.subject_to( ... );
  
end

% ---- path constraints -----------

limit = track.maxspeed;
opti.subject_to(speed  <=   limit(pos) + epsilon_speed); % track speed limit
opti.subject_to(0 <= U <= 1);  % control is limited

% ---- boundary conditions --------
opti.subject_to(pos(1)==pos0);   % use initial position
opti.subject_to(speed(1)==v0); % use initial speed
opti.subject_to(epsilon_speed >= 0); % slack lower bound

% Pass parameter values
opti.set_value(pos0, 0.0);
opti.set_value(v0, 0);

% ---- Setup solver NLP    ------
ops = struct;
ops.ipopt.print_level = 0;
ops.ipopt.tol = 1e-3;
opti.solver('ipopt', ops);
sol = opti.solve();   % actual solve

% ---- Plot predicted trajectory for debugging ------

figure(1); clf
title('Predicted trajectory from initial state')
hold on; grid on
plot(limit(sol.value(pos)),'k--','linewidth',2);
plot(sol.value(pos),'linewidth',2);
plot(sol.value(speed),'linewidth',2);
stairs(1:N, sol.value(U(1,:)),'linewidth',2);
stairs(1:N,-sol.value(U(2,:)),'linewidth',2);
legend('speed limit','pos','speed','throttle','brake','Location','northwest')
