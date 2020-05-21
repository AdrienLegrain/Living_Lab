function Q=HeatELG(T_ext)
%HEATELL Summary of this function goes here
%   Detailed explanation goes here
p1 = -1.664143004704622e+03;
p2 = 3.091314730452895e+04;

if(abs(p1*T_ext)>=abs(p2))
    Q=0;
else
    Q = p1*T_ext+p2;
end

end