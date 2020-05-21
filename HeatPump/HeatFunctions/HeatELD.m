function Q=HeatELD(T_ext)
%HEATELL Summary of this function goes here
%   Detailed explanation goes here
p1 = -2.645835944294058e+03;
p2 = 5.440957052252987e+04;
if(abs(p1*T_ext)>=abs(p2))
    Q=0;
else
    Q = p1*T_ext+p2;
end

end