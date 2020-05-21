function Q=HeatELE(T_ext)
%HEATELL Summary of this function goes here
%   Detailed explanation goes here
p1 = -2.618620840169914e+03;
p2 = 5.074286506169281e+04;
if(abs(p1*T_ext)>=abs(p2))
    Q=0;
else
    Q = p1*T_ext+p2;
end

end