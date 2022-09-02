function dydt = Function(t,y)
% For use in ODE45. Not much here.

global Eff Gee 
dydt = [Eff(y(1),y(2)); Gee(y(1),y(2))];
end

