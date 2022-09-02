function [compactpropensity] = compact_propensity(monomials,reactionrates,s1,s2,s3,s4)
% Takes an array of monomials and reactionrates, returns a symbolic
% 2-dimensional vector with the correct form of propensity for modified
% AutoGillespie. 

syms skeletonpropensity [1 length(monomials)]
for i = 1 : length(monomials)
    mono = monomials(i);
    
    if polynomialDegree(mono) == 0
    result = mono;
    else 
        s1propensity = Gillespieloop(s1,polynomialDegree(mono,s1));
        s2propensity = Gillespieloop(s2,polynomialDegree(mono,s2));
        s3propensity = Gillespieloop(s3,polynomialDegree(mono,s3));
        s4propensity = Gillespieloop(s4,polynomialDegree(mono,s4));
        result = s1propensity*s2propensity*s3propensity*s4propensity;
    end
    skeletonpropensity(i) = result;
end

compactpropensity = [sum(skeletonpropensity .* eval(reactionrates>=0) .* abs(reactionrates)),...
    sum(skeletonpropensity .* eval(reactionrates<0).* abs(reactionrates))];
end

