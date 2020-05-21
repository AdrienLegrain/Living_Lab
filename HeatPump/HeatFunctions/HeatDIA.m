function Q = HeatDIA(T_ext)
%HEATDIA Summary of this function goes here
%   Detailed explanation goes here
p1 = -9.808911274826511e+02;
p2 = 1.816758504466745e+04; 
if(abs(p1*T_ext)>=abs(p2))
    Q=0;
else
    Q = p1*T_ext+p2;
end
end

