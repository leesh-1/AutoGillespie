function [x1,x2,x3,x4] = Modified_Gillespie(Propensity_Handle,t,initial_condition)
% Impliments Gillespie with 8 reaction channels

% Initialize 
x1 = zeros(size(t));
x2 = zeros(size(t));
x3 = zeros(size(t));
x4 = zeros(size(t));

% Total Reaction Time
T = 0;

% Initial Condition, compute properly
x1(1) = initial_condition(1);
x2(1) = initial_condition(2);
x3(1) = initial_condition(3);
x4(1) = initial_condition(4);
j = 2; x1(j)=x1(j-1); x2(j)=x2(j-1); x3(j)=x3(j-1); x4(j)=x4(j-1);

while T < t(j)
    r1 = rand; r2 = rand;

    % Propensity computation
    prop = Propensity_Handle(x1(j),x2(j),x3(j),x4(j));
    % Next reaction time
    tau = 1/sum(prop) * log(1/r1);
    T = T + tau;
    while j < length(t) && t(j) <= T
        j = j+1;
        x1(j) = x1(j-1);
        x2(j) = x2(j-1);
        x3(j) = x3(j-1);
        x4(j) = x4(j-1);
    end
    if j > length(t) || T > t(end)
        break
    end

    % Update reaction effects
    prop1 = prop./sum(prop);
    if r2 < sum(prop1(1:1))
        x1(j) = x1(j)+1;
    elseif r2 < sum(prop1(1:2))
        x1(j) = x1(j)-1;
    elseif r2 < sum(prop1(1:3))
        x2(j) = x2(j)+1;
    elseif r2 < sum(prop1(1:4))
        x2(j) = x2(j)-1;
    elseif r2 < sum(prop1(1:5))
        x3(j) = x3(j)+1;
    elseif r2 < sum(prop1(1:6))
        x3(j) = x3(j)-1;
    elseif r2 < sum(prop1(1:7))
        x4(j) = x4(j)+1;
    else
        x4(j) = x4(j)-1;
    end
end
end