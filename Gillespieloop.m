function [propensityelement] = Gillespieloop(s,monomialDegree)
% To properly compute Propensity functions.

propensityelement = s;
if monomialDegree == 0
    propensityelement = 1;
else
    for deg = 1 : monomialDegree-1
        propensityelement = propensityelement * (s-deg);
    end
end
end

