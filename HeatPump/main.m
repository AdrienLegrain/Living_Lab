clear all
close all

%% MED Data Plot
ImportData
Time = heatdata{:,2}
Heat = heatdata{:,3}

figure 
plot(Time,Heat,'.')
xlabel('Time')
ylabel('Heat demand of MEB')

%% Heat pump Model
% Parameters
K = 1 %Efficiency of the compressor

% T_cin = 1 % Temp.  Of refrig.  Into Condenser (°C)
% T_cout = 1 %Temp.  of refrig.  out of Condenser (°C)
% T_ein = 1 % Temp.  of refrig.  into Evaporator (°C)
% T_eout = 1% Temp.  of refrig.  out of Evaporator (°C)
% 
% T_cin = 47
T_ein= 7;
dm_w1 = 11.38; % Massflow rate of refrigerant in condenser (kg/s)
C_e = 4180 ; % Heat capacity of the evaporator (J/kg K)
m_e = 1 ; % mass of evaporator
dm_w2 = 11.38 ; % Massflow rate of refrigerant in condenser (kg/s)
C_c = 4180 ; % Heat capacity of the Condenser (J/kg K)
m_c = 1 ; % mass of condenser
C_f = 1.107 ; % Heat capacity of the refrigerant (J/kg K)
dW = 5 ; % dif of the work done to the HP cycle -> Energy input in Joules
tspan = [0 10];
m_water = dm_w2*tspan(2);
Cwater = 4180;
CairVbuilding = 10000;
% Equations 

Ke = dm_w1*C_f/C_e/m_e;
Kc = dm_w2*C_f/C_c/m_c;
Hd = 411.910431315784;
%for systems of ODE see: https://ch.mathworks.com/help/symbolic/solve-a-system-of-differential-equations.html
% x = [Teout,Teout',Tcout, Tcout',Tcin,Tcin'] 
% eq1 = -Ke*x(1)-Ke*x(5) - dW/(C_e*m_e)*((x(2)+2*273.15)/(x(2)-T_ein-x(1))-1)
% eq2 = -Kc*x(2)+Kc*x(5) - dW/(C_c*m_c)*((x(2)+2*273.15)/(x(2)-T_ein-x(1))-1)
% eq3 = (Cwater*m_w2*(x(3)-x(5))-Hd)/(Cair*Vbuilding)
% x(5) = (x(3)*m_water*Cwater-Hd)/(CairVbuilding + m_water*Cwater)
f = @(t,x) [x(2);-Ke*x(1)+Ke*x(5) - dW/(C_e*m_e)*((x(2)+2*273.15)/(x(2)-T_ein-x(1))-1); x(4);-Kc*x(2)+ Kc*x(5) - dW/(C_c*m_c)*((x(2)+2*273.15)/(x(2)-T_ein-x(1))-1);x(5);(x(3)*m_water*Cwater-Hd)/(CairVbuilding + m_water*Cwater)];
tspan = [0 10];
X0 = [47; 0; 35; 0; 47; 0] 
[T,X] = ode45(f,tspan,X0)
Sol = ode45(f,tspan,X0)


%%
% for dW=1:1
%     sol = Heatpump(K,T_cin, T_ein, dm_w1, dm_w2, C_e, C_c, m_e, m_c, C_f, dW, tspan, X0)
% 
%     test1 = linspace(0,tspan(2),250);
%     test2 = deval(sol,test1);
% 
%     figure (1)
%     plot(test1,test2(1,:))
%     hold on
% end
% %%
% function [Tcin] = THeatout(m_water,Hd,Tcout,Teout)
% Cwater = 4180
% CairVbuilding = 10000
% Tcin = (Tcout*m_water*Cwater-Hd)/(CairVbuilding + m_water*Cwater)
% end
% m_water = dm_w2*4
% Tcin = THeatout(1000,411,X(1,3), X(1,1))
% 
