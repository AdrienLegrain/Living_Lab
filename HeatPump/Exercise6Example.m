clear; close all; clc;

car_dynamics

%% Part 1: Define and test your integrators
addpath('./integrators')

%%% Part 1
%%% Define your integrator functions in the ./integrators folder

%%% Part 1 - 1 Defines f_{discrete} here
h = 0.01; % integration time step
f_discrete_RK4 = @(x,u) RK4(x,u,h,f);
f_discrete_euler = @(x,u) Euler(x,u,h,f);

%%% Part 1-2 Test the integration call
X0 = [0;0.5]; U0 = [0;-0.1];
try
  Xh_RK4 = f_discrete_RK4(X0,U0)
  Xh_euler = f_discrete_euler(X0,U0)
catch
  warning('You need to implement the functions RK4 and Euler in the integrators directory');
end

%%% If everything is correct, output should be
%%% Xh_RK4 = [0.0050,0.5033] and Xh_euler = [0.0050,0.5033]
clear; close all; clc;

car_dynamics

%% Part 2: Simulate the system

h = 1; % Set sample period

t = 0:h:10; % Sample times
X0 = [0;0.5]; % Initial state
U0 = 0; % Initial input
Uref = [0.5+0.5*sin(t);zeros(size(t))]; % Sample the input function at the sample period

%%% Simulate using ODE45
ode.name = 'ODE45';
ode.f_discrete = @(X,U) ode45(@(t,x) f(x,U),[0 h], X');
ode.X = X0;
for k=1:length(t)-1
  res = ode.f_discrete(ode.X(:,k),Uref(:,k));
  ode.X(:,k+1) = res.y(:,end);
end

%%% Simulate using RK4

rk4.name = 'RK4';
% --->>> Code here
rk4.f_discrete = @(x,u) RK4(x,u,h,f);
rk4.X = X0;
for k = 1:length(t)-1
  rk4.X(:,k+1) = rk4.f_discrete(rk4.X(:,k), Uref(:,k));
end
% --->>> Code here

%%% Simulate using Euler

eur.name = 'Euler';
% --->>> Code here
eur.f_discrete = @(x,u) Euler(x,u,h,f);
eur.X = X0;
for k = 1:length(t)-1
  eur.X(:,k+1) = eur.f_discrete(eur.X(:,k), Uref(:,k));
end
% --->>> Code here


%% Plot results
plotSims(t, {ode, rk4, eur});


clear; close all; clc;

car_dynamics
h = 0.5;

% Load your RK4 integrator from Ex 1
addpath('./integrators')
f_discrete = @(x,u) RK4(x,u,h,f);

%% Exercise 2.1: Linearizing a nonlinear function

%%% Compute the jacobians using algorithmic differentiation
% Note syntax for using casadi from this template for future use
import casadi.*   % import casadi once before using

X0 = SX.sym('X0',2,1); % declare a symbolic variable X0 of size 2x1
U0 = SX.sym('U0',2,1); % declare a symbolic variable U0 of size 2x1

x_next = f_discrete(X0,U0); % Calculate the next step symbolically

A_jac = jacobian(x_next,X0); % returns jacobian for an expression (x_next) with respect to X0
B_jac = jacobian(x_next,U0); % returns jacobian for an expression (x_dot) with respect to U0

% convert A_algorithmic, B_algorithmic expressions to a callable functions
A_algorithmic = casadi.Function('A_algorithmic',{X0,U0},{A_jac});
B_algorithmic = casadi.Function('A_algorithmic',{X0,U0},{B_jac});

%% Compare the algorithmic results to finite differences

X0 = [1;0.5]; U0 = [0;0]; % assign some numeric values now to compute jacobians using the function

fprintf('A_algorithmic(X0,U0)=\n')
disp(full(A_algorithmic(X0,U0)))
fprintf('B_algorithmic(X0,U0)=\n')
disp(full(B_algorithmic(X0,U0)))

try
  fprintf('A_finite_diff(X0,U0)=\n')
  disp(jac_x(X0,U0,f_discrete))
  fprintf('B_finite_diff(X0,U0)=\n')
  disp(jac_u(X0,U0,f_discrete))
catch
  warning('Implement the functions jac_x and jac_u to compute the jacobians using finite differences');
end
clear; close all; clc;

Ex7_2_1 % Load your linearization routines

%%% ---> Your code here
% Maps from X,U => next state
f_linear_algo = @(X,U) error('Implement me');  % Linearized system using casadi
f_linear_fd   = @(X,U) error('Implement me');  % Linearized system using finite differences

A_al = full(A_algorithmic(X0,U0));
B_al = full(B_algorithmic(X0,U0));
f_linear_algo = @(X,U) A_al * (X - X0) + B_al * (U - U0) + f_discrete(X0,U0);  % Linearized system using casadi

A_fd = jac_x(X0,U0,f_discrete);
B_fd = jac_u(X0,U0,f_discrete);
f_linear_fd   = @(X,U) A_fd * (X - X0) + B_fd * (U - U0) + f_discrete(X0,U0);  % Linearized system using finite differences

%%% ---> Your code here


%% Simulate

sim_algo.name = 'Linearized with algorithmic differentiation';
sim_algo.f = f_linear_algo;

sim_fd.name = 'Linearized with finite differences';
sim_fd.f = f_linear_fd;

sim_rk4.name = 'Nonlinear RK4';
sim_rk4.f = f_discrete;

sims = {sim_algo, sim_fd, sim_rk4};

t = 0:h:10; % Sample times
X0 = [0;0.5]; % Initial state
U0 = 0; % Initial input
Uref = [0.5+0.5*sin(t);zeros(size(t))]; % Sample the input function at the sample period

for i = 1:length(sims)
  try
    sims{i}.X = X0;
    for k = 1:length(t)-1
      sims{i}.X(:,k+1) = sims{i}.f(sims{i}.X(:,k), Uref(:,k));
    end
  catch
    warning('Implement the function for %s', sims{i}.name);
  end
end

%% Plot results
plotSims(t, sims);

%% Compute linearized dynamics

A_fd = jac_x(X0,U0,f);
B_fd = jac_u(X0,U0,f);

A = full(A_auto_function(X0,U0));
B = full(B_auto_function(X0,U0));

f_lin_fd   = @(X,U) ... ; % using finite difference
f_lin_auto = @(X,U) ... ; % using algorithmic derivatives

%% Simulate
h1 = 0.1;
Uref = @(t) [0.5+0.5*sin(t);0];

N_h1 = floor(10/h1);
f_discrete_auto = @(x,u) RK4(x,u,h1,f_lin_auto);
f_discrete_fd = @(x,u) RK4(x,u,h1,f_lin_fd);
f_discrete = @(x,u) RK4(x,u,h1,f);

X = []; X_auto = []; X_fd = []; t=[];
% Run for h1 
x0 = X0; x0_auto = X0; x0_fd = X0; t0 = 0; % initialize
for k=1:N_h1-1
% simulate the automatic differentiation system
Xstep = f_discrete_auto(x0_auto,Uref(t0));
X_auto=[X_auto,Xstep]; 
% simulate the finite difference system
Xstep = f_discrete_fd(x0_fd,Uref(t0));
X_fd=[X_fd,Xstep];
% simulate the nonlinear system
Xstep = f_discrete(x0,Uref(t0));
X=[X,Xstep]; 
t = [t,t0+h1]; % collect results for plotting
% update for next iteration
x0_auto = X_auto(:,end); 
x0_fd = X_fd(:,end);
x0 = X(:,end);
t0 = t(end); 
end

% plot the trajectories
figure(2)
subplot(1,2,1)
l0 = plot(t,X(1,:),'linewidth',2);
xlabel('t (s)')
ylabel('\lambda')
subplot(1,2,2)
plot(t,X(2,:),'linewidth',2);
xlabel('t (s)')
ylabel('v')

% plot the RK4 results


subplot(1,2,2)
hold on
plot(t,X_auto(2,:),'linewidth',2);
figure(2)
subplot(1,2,2)
hold on
plot(t,X_fd(2,:),'--','linewidth',2);

figure(2)
subplot(1,2,1)
hold on
l21 = plot(t,X_auto(1,:),'linewidth',2);

subplot(1,2,1)
hold on
l22 = plot(t,X_fd(1,:),'--','linewidth',2);
legend([l0,l21,l22],'nonlinear','auto. diff.','finite diff.','location','northwest')



    clear; close all;
import casadi.*
addpath('./plottingcode')

car_dynamics

%%%%% REMOVE
% h = 0.1;
h = 0.025;
f_discrete = @(x,u) RK4(x,u,h,f);
%%%%% REMOVE

car_dynamics

%%% Choose your track
track = generate_track('simple');
% track = generate_track('complex');

%% Ex 3 MPC problem setup

opti = casadi.Opti(); % Optimization problem

N = 10; % MPC horizon

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

% % Minimize braking
% opti.minimize(...
%   -10*sum(X(2,:))  + ... % Max velocity
%   0.1*U(1,:)*U(1,:)' + ... % Minimize accel
%   1000*U(2,:)*U(2,:)'   + ... % Minimize braking
%   10000*(epsilon_speed(1,:)*epsilon_speed(1,:)' + sum(epsilon_speed))); % Soft constraints


% ---- multiple shooting --------
for k=1:N % loop over control intervals
  
  %%%%% REMOVE
  opti.subject_to(X(:,k+1) == f_discrete(X(:,k), U(:,k)));
  %%%%% REMOVE
  
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

%% closed loop simulation

sim.x = [0;0]; % Start at position zero with zero speed
sim.u = [];
sim.t = 0;

figure(2); clf; plot_track(track); hndl = [];

i = 0;
while sim.x(1,end) < 1  %% Take one loop around the track
  i = i + 1;
  sim.t(end+1) = sim.t(end) + h; 
  
  % set position and speed variables to their current values
  opti.set_value(X0, sim.x(:,i));
  
  % ---- call the MPC controller ------
  sol = opti.solve();
  u_mpc = sol.value(U(:,1)); % Get the MPC input
  
  % ---- Simulate the system
  next = ode45(@(t,x) f(x,u_mpc), [0, h], sim.x(:,i));
  sim.x(:,i+1) = next.y(:,end);
  sim.u(:,i) = u_mpc;
  
  % ---- Plot the car (delete this line if it's taking too long)
  hndl = plot_car(track, sol.value(pos(1)), sol.value(pos(2:end)), hndl);
  %   pause(0.05); drawnow
  
  % warm start the next iteration
  opti.set_initial(U, sol.value(U));
  opti.set_initial(X, sol.value(X));
  opti.set_initial(epsilon_speed, sol.value(epsilon_speed));
end

%% Plots for the closed loop results
figure(3); clf; hold on; grid on;
plot(sim.t, limit(sim.x(1,:)),'k--', 'linewidth',4); % Speed constraint
plot(sim.t, sim.x(1,:), 'linewidth',2); % Position
plot(sim.t, sim.x(2,:), 'linewidth',2); % Speed

stairs(sim.t(1:end-1), sim.u(1,:),'linewidth',2); % Acceleration
stairs(sim.t(1:end-1),-sim.u(2,:),'linewidth',2); % Braking

legend('speed limit','pos','speed','throttle','brake','Location','northwest')
title('Closed loop trajectory')

xlabel('Time (s)')
