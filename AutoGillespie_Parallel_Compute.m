function [xx_1,xx_2,xx_3,xx_4] = AutoGillespie_Parallel_Compute(polyfunc1,polyfunc2,Gilltime,initial_condition,NumRealize)
% Impliments the AutoGillespie algorithm, with parallel computing enabled.

% Initialize
xx_1 = zeros(length(Gilltime),NumRealize);
xx_2 = zeros(length(Gilltime),NumRealize);
xx_3 = zeros(length(Gilltime),NumRealize);
xx_4 = zeros(length(Gilltime),NumRealize);

% Apply QSST and return the equations in vectorized format
[monomials1,reactionrates1,monomials2,reactionrates2,monomials3,...
    reactionrates3,monomials4,reactionrates4,s1,s2,s3,s4]...
    = Impliment_QSST(polyfunc1,polyfunc2);

% Derive a propensity vector, which will be used for Gillespie Simulation
[Propensity_Handle] = Return_Compact_Propensity(monomials1,reactionrates1,...
    monomials2,reactionrates2,monomials3,reactionrates3,monomials4,reactionrates4...
    ,s1,s2,s3,s4,Gilltime,initial_condition);

% Run the (modified) Gillespie algorithm with compact propensity functions
parfor i = 1 : NumRealize
    [xx_1(:,i),xx_2(:,i),xx_3(:,i),xx_4(:,i)] = Modified_Gillespie(Propensity_Handle,Gilltime,initial_condition);
end
end