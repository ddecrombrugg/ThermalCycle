function [Y] = adjust(T,h,s)
% The aim of this function is to complete the XSteam function. In fact, we cannot calculate
%       - the enthalpy in function of the temperature and the entropy
%       - the entropy in function of the temperature and the enthalpy
% INPUTS
%       - T : Temperature of the research state [°C]
%       - h : Enthalpy of the research state [kJ/kg]
%           (vector of 2 values if it is the researched value)
%       - s : Entropy of the research state [kJ/(kg*K)]
%           (vector of 2 values if it is the researched value)
% OUTPUT
%       - Y : Whether the enthalpy or the entropy of the researched value
%

tol = 1e-3; nmax = 5e2;

if h == 0
    [Y1,Y2] = temp_limits_s(s);
elseif s == 0
    [Y1,Y2] = temp_limits_h(h);
else
    fprintf('Whether vector ''h'' or vector ''s'' must be unknown.\n')
    Y = NaN; return;
end


err = 1; n = 0;
while (abs(err) > tol && n < nmax)
    Y = (Y1 + Y2)/2;
    if h == 0
        err = T - XSteam('T_hs',Y,s);
        if err > 0
            Y1 = Y;
        else
            Y2 = Y;
        end
    else
        err = T - XSteam('T_hs',h,Y);
        if err < 0
            Y1 = Y;
        else
            Y2 = Y;
        end
    end
    n = n +1;
end

end


function [h1,h2] = temp_limits_s(s)
    h = 0:5e3;
    T = zeros(length(h),1);
    for i = 1:length(h)
        T(i) = XSteam('T_hs',h(i),s);
    end
    h1 = h(find(T == min(T)));
    h2 = h(find(T == max(T)));
end

function [s1,s2] = temp_limits_h(h)
    s = 0:1e-2:10;
    T = zeros(length(s),1);
    for i = 1:length(s)
        T(i) = XSteam('T_hs',h,s(i));
    end
    s1 = s(find(T == max(T)));
    s2 = s(find(T == min(T)));
end

