function dydt = FunctionTaylor(t,y)
% For use in ODE45. Not much here.

global EffTay GeeTay 
dydt = [EffTay(y(1),y(2)); GeeTay(y(1),y(2))];
end

