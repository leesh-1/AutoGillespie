function [Propensity_Handle] = Return_Compact_Propensity(monomials1,reactionrates1,...
    monomials2,reactionrates2,monomials3,reactionrates3,monomials4,reactionrates4...
    ,s1,s2,s3,s4,t,initial_condition)
% Returns "Compact Propensity" functions which enable the formalization of
% any 2-dimensional polynomial of any degree into a pseudo-reaction network
% of just 8 reaction channels, under the canonical inversion.
propensity1 = compact_propensity(monomials1,reactionrates1,s1,s2,s3,s4);
propensity2 = compact_propensity(monomials2,reactionrates2,s1,s2,s3,s4);
propensity3 = compact_propensity(monomials3,reactionrates3,s1,s2,s3,s4);
propensity4 = compact_propensity(monomials4,reactionrates4,s1,s2,s3,s4);
FullPropensity = [propensity1,propensity2,propensity3,propensity4];
Propensity_Handle = matlabFunction(FullPropensity,'Vars',{s1 s2 s3 s4});
end

