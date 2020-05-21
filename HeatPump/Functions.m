
% function Q = Heat
%   Heat.ELL=@ELL;
% %   Heat.ELA=@ELA;
% end

% function Q= ELL(T_ext)
%       Q = FitsValue(1,1)*T_ext+FitsValue(1,2);
% end

ELL = @(T_ext) (FitsValue(1,1)*T_ext+FitsValue(1,2));

% function Q= ELA(T_ext)
%       Q = FitsValue(2,1)*T_ext+FitsValue(2,2);
% end

% 