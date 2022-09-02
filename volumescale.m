function [result] = volumescale(V,degree)
% Volume scaling, for use by Gillespie algorithm.
if degree == 0
    result = V;
else
    result = V^(-degree+1);
end
end

