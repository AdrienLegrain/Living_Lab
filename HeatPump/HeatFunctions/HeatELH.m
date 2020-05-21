function Q=HeatELH(T_ext)
%HEATELL Summary of this function goes here
%   Detailed explanation goes here
p1 = -9.484072345530843e+02;
p2 = 1.853116597446906e+04;
if(abs(p1*T_ext)>=abs(p2))
    Q=0;
else
    Q = p1*T_ext+p2;
end

end