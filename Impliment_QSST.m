function [monomials1,reactionrates1,monomials2,reactionrates2,monomials3,...
    reactionrates3,monomials4,reactionrates4,s1,s2,s3,s4] = Impliment_QSST(polyfunc1,polyfunc2)
% Applies Quasi-Steady State Transformation to polyfunc1, polyfunc2,
% returns the general system in a vectorized format.

% Required parameters for Quasi-Steady State Transformation
global omega p1 p2 mu

% This is a technicality in the code. The individual functions composing
% this code is capable of forming another automated simulation algorithm,
% in which case volume scaling is important. For this file, we need not
% scale, so let volume V = 1. That is, we have already scaled the ODE
% function in the main file Hopf, so we do not scale twice.
V=1;

syms s3 s4 ghost kg ghost1 kg1

% Extract the two symbolic variables in polyfunc1, polyfunc2
f = abs(polyfunc1) + abs(polyfunc2);
if isempty(f)
    polyfunc1_and_polyfunc2_havenosymbolicvariables_error = 1
end
s = symvar(f,2);
if isempty(s)
    syms ghost2 ghost3
    s = [s,ghost2,ghost3];
elseif length(s) == 1
    syms ghost2
    s = [s,ghost2];
end
s1 = s(1); s2 = s(2);

% Add in ghost variables to ensure proper implimentation
polyfunc1 = polyfunc1 + kg * ghost + kg1 * ghost1;
polyfunc2 = polyfunc2 + kg * ghost + kg1 * ghost1;
polyfunc1 = expand(polyfunc1);
polyfunc2 = expand(polyfunc2);

% Split up into polyfunc1, polyfunc2 into monomials and reaction rates
[monomials1,reactionrates1] = Gillespie_format(polyfunc1,s(1),s(2),s3,s4,ghost,ghost1,V);
[monomials2,reactionrates2] = Gillespie_format(polyfunc2,s(1),s(2),s3,s4,ghost,ghost1,V);

% Knock out Ghost Variables
kg = 0; kg1 = 0;
SubbedReacRate1 = subs(reactionrates1);
SubbedReacRate2 = subs(reactionrates2);
% Detect if polyfunc1, polyfunc2 have cross-negative terms
crossneg_tester1 = matlabFunction(monomials1 .* eval(SubbedReacRate1 < 0),'Vars',{s(1) s(2)});
crossneg_tester2 = matlabFunction(monomials2 .* eval(SubbedReacRate2 < 0),'Vars',{s(1) s(2)});

% Obtain Pointers to cross-negative terms
replaceposition1 = crossneg_tester1(0,1) ~= zeros(size(SubbedReacRate1));
replaceposition2 = crossneg_tester2(1,0) ~= zeros(size(SubbedReacRate2));

% We replace the monomials specified in the above positions with QSST
newmonomials1 = monomials1.*replaceposition1*s(1)*s3*omega^(-1)*p1;
monomials1 = (replaceposition1 == 0).* monomials1 + newmonomials1;
newmonomials2 = monomials2.*replaceposition2*s(2)*s4*omega^(-1)*p2;
monomials2 = (replaceposition2 == 0).* monomials2 + newmonomials2;

% The general system is now kinetic. We represent general system first.
% Ghosts must be added in once more
syms kg kg1

% In case these are not assigned during inversion process, i.e. no
% cross-negative terms are detected
monomials3 = sym(0); reactionrates3 = sym(0);
monomials4 = sym(0); reactionrates4 = sym(0);

if ~isequal(replaceposition1,zeros(size(replaceposition1)))
    polyfunc3 = 1/mu * (omega - s(1)*s3*p1);
    polyfunc3 = polyfunc3 + kg * ghost + kg1 * ghost1;
    polyfunc3 = expand(polyfunc3);
    [monomials3,reactionrates3] = Gillespie_format(polyfunc3,s(1),s(2),s3,s4,ghost,ghost1,V);
end

if ~isequal(replaceposition2,zeros(size(replaceposition2)))
    polyfunc4 = 1/mu * (omega - s(2)*s4*p2);
    polyfunc4 = polyfunc4 + kg * ghost + kg1 * ghost1;
    polyfunc4 = expand(polyfunc4);
    [monomials4,reactionrates4] = Gillespie_format(polyfunc4,s(1),s(2),s3,s4,ghost,ghost1,V);
end

% We have vectorized the general system, so now we feed it into Gillespie.
% We knock out the remaining ghost variables
kg = 0; kg1 = 0; ghost = 0; ghost1 = 0;
reactionrates1 = subs(reactionrates1);
reactionrates2 = subs(reactionrates2);
monomials1 = subs(monomials1);
monomials2 = subs(monomials2);
reactionrates3 = subs(reactionrates3);
reactionrates4 = subs(reactionrates4);
monomials3 = subs(monomials3);
monomials4 = subs(monomials4);
end
