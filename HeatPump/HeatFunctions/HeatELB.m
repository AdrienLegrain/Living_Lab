function Q=HeatELB(T_ext)
%HEATELL Summary of this function goes here
%   Detailed explanation goes here
p1 = -2.195539601430054e+03;
p2 = 4.292950695565090e+04;
if(abs(p1*T_ext)>=abs(p2))
    Q=0;
else
    Q = p1*T_ext+p2;
end
end
