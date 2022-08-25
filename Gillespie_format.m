function [monomials,scaledreactionrates] = Gillespie_format(polyfunc,s1,s2,s3,s4,ghost,ghost1,V)
% Takes EXPANDED polynomial function with NO BRACKETS, reshapes into 
% Gillespie form. Functions need to be expanded to use "children". 

% s(1),s(2),s3,s4,ghost,ghost1 are all sybmolic variables.
NumTerms = numel(children(polyfunc));

syms scaledreactionrates [1 NumTerms]
syms reactionrates [1 NumTerms]
syms monomials [1 NumTerms]
for j = 1 : NumTerms
    monomials(j) = children(polyfunc,j);
    reactionrates(j) = coeffs(monomials(j), [s1,s2,s3,s4,ghost,ghost1]);
    scaledreactionrates(j) = reactionrates(j)*volumescale(V,polynomialDegree(monomials(j)));
end
monomials = monomials ./ reactionrates;
end

