function [Sol] = Heatpump(K,T_cin, T_ein, dm_w1, dm_w2, C_e, C_c, m_e, m_c, C_f, dW, tspan, X0)
%Heatpump 
%   Detailed explanation goes here
[K,T_cin, T_ein, dm_w1, dm_w2, C_e, C_c, m_e, m_c, C_f]

Ke = dm_w1*C_f/C_e/m_e
Kc = dm_w2*C_f/C_c/m_c

f = @(t,x) [x(2);-Ke*x(1)+Ke*T_cin - dW/(C_e*m_e)*((x(2)+2*273.15)/(x(2)-T_ein-x(1))-1); x(4);-Kc*x(2)+ Kc*T_cin - dW/(C_c*m_c)*((x(2)+2*273.15)/(x(2)-T_ein-x(1))-1)];
Sol = ode45(f,tspan,X0)
end

