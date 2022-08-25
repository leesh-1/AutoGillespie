function [x,yuse,u,v] = get_vectorgradients(normalize,xnum,ynum,xlow,xhigh,ylow,yhigh,f)
% Obtains vector gradients of function f, at points specified by y1, y2
% below. If normalize = 1 in input, the function normalizes the vector
% gradients.

% Initialization and set-up
y1 = linspace(xlow,xhigh,xnum);
y2 = linspace(ylow,yhigh,ynum);
[x,yuse] = meshgrid(y1,y2);
u = zeros(size(x));
v = zeros(size(x));

% We want the gradients at the start time t=0
t=0;
for i = 1:numel(x)
    Yprime = f(t,[x(i); yuse(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end
if normalize == 1
    for i = 1:numel(x)
        Vmod = sqrt(u(i)^2 + v(i)^2);
        u(i) = u(i)/Vmod;
        v(i) = v(i)/Vmod;
    end
end
end