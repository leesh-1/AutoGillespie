% A sketch of the AutoGillespie Algorithm. A sample application to a
% bifurcated system demonstrating a Hopf Bifurcation. The code has been
% altered significantly for clarity, and therfore may no longer be suitable
% for general use. This version only operates on 2D systems. The code needs
% to be modified when there are more or less than 2 symbolic variables in 
% the ODE system.

clear
close all
clc


global V Eff Gee EffTay GeeTay translatebyx translatebyy omega p1 p2 mu

% Volume of reacting solution, for the chemical simulation
V = 500;

% We make an affine transformation (a diagnal push) to relocate the
% dynamics of interest to the positive quadrant.
translatebyx = 6;
translatebyy = 6;

% We Taylor expand to this order, the RHS of the input ODE system.
expand_to_order = 4;

% For better illustration, we will plot the phase plane. Set axis bounds.
width = translatebyx;
xlow = (-width+translatebyx)*V;
xhigh = (width+translatebyx)*V;
ylow = (-width+translatebyy)*V;
yhigh = (width+translatebyy)*V;

figure(1)
axis([xlow xhigh ylow yhigh])
hold on
title('Original Phase Plane','Fontsize',16)
xlabel('Specimen X Molecule Count','Fontsize',16)
ylabel('Specimen Y Molecule Count','Fontsize',16)

figure(2)
axis([xlow xhigh ylow yhigh])
hold on
title('Taylor Induced Phase Plane','Fontsize',16)
xlabel('Specimen X Molecule Count','Fontsize',16)
ylabel('Specimen Y Molecule Count','Fontsize',16)

% We trace deterministic dynamics in the phase plane. We divide the
% vertical/horizontal axis into ynum/xnum data points.
xnum = 15;
ynum = 15;

% Initial conditions for deterministic plots on phase plane
y1axis = linspace(xlow,xhigh,xnum);
y2axis = linspace(ylow,yhigh,ynum);
[Y1,Y2] = meshgrid(y1axis,y2axis);


% We give the ODE system we wish to chemically simulate
syms x1 y1
% Give parameters
xi = 10; zeta = 10;
% RHS of 2D ODE system we wish to chemically simulate
Eff = xi*x1 - zeta*y1 - x1*(x1^2+y1^2)/V^2;
Gee = zeta*x1 + xi*y1 - y1*(x1^2+y1^2)/V^2;

% We Taylor expand the system at this point
pnt_of_expand = [0 0];
Taylorf = taylor(Eff,[x1 y1],pnt_of_expand,'Order',expand_to_order);
Taylorg = taylor(Gee,[x1 y1],pnt_of_expand,'Order',expand_to_order);

% We make the affine transformation to push dynamics to positive cone
Eff = subs(Eff,[x1,y1],[x1-translatebyx*V,y1-translatebyy*V]);
Gee = subs(Gee,[x1,y1],[x1-translatebyx*V,y1-translatebyy*V]);
Taylorf = subs(Taylorf,[x1,y1],[x1-translatebyx*V,y1-translatebyy*V]);
Taylorg = subs(Taylorg,[x1,y1],[x1-translatebyx*V,y1-translatebyy*V]);

% Return as a function handle
EffTay = matlabFunction(Taylorf,'Vars',{x1 y1});
GeeTay = matlabFunction(Taylorg,'Vars',{x1 y1});
Eff = matlabFunction(Eff,'Vars',{x1 y1});
Gee= matlabFunction(Gee,'Vars',{x1 y1});

% Draw Phase Planes
tstart = 0; tend = 20;
for n = 1 : numel(Y1)

    [t,y] = ode45(@Function,[tstart tend],[Y1(n); Y2(n)]);
    figure(1)
    plot(y(:,1),y(:,2),'c')

    [t_0,y_0] = ode45(@FunctionTaylor,[tstart tend],[Y1(n); Y2(n)]);
    figure(2)
    plot(y_0(:,1),y_0(:,2),'c')

end

%% Plot streamlines for Phase Plane

% The more data points, the better the results are
xnum = 120; ynum = 120;

% See get_vectorgradients function; we do NOT normalize the gradients
normalize = 0;

f = @Function;
[x,yuse,u,v] = get_vectorgradients(normalize,xnum,ynum,xlow,xhigh,ylow,yhigh,f);
figure()
hold on
streamslice(x,yuse,u,v)
axis image
title('Original Phase Plane','Fontsize',16)
xlabel('Specimen X Molecule Count','Fontsize',16)
ylabel('Specimen Y Molecule Count','Fontsize',16)

g = @FunctionTaylor;
[x,yuse,u,v] = get_vectorgradients(normalize,xnum,ynum,xlow,xhigh,ylow,yhigh,g);
figure()
hold on
streamslice(x,yuse,u,v)
axis image
title('Taylor Induced Phase Plane','Fontsize',16)
xlabel('Specimen X Molecule Count','Fontsize',16)
ylabel('Specimen Y Molecule Count','Fontsize',16)

%% Now use AutoGillespie for Chemical Simulation

% Define the timestepping
Gilltime = linspace(0,2,10^7);
% Average over this many trajectories
NumRealize = 1;

% Inputs for AutoGillespie
initial_condition = [2256 1645 1 1];
polyfunc1 = Taylorf;
polyfunc2 = Taylorg;

% Deterministic solution to compare with stochastic solution
[t_9,y_9] = ode45(@Function,Gilltime,[initial_condition(1); initial_condition(2)]);

% For Quasi-Steady State Transformation
omega = 1; % omega_x1 = omega_y1 = 1
p1 = 1; p2 = 1; % p_x1(x1,y1) = p_x2(x1,y1) = 1
mu = 10^(-6);

% Using the notation in the dissertation, xx_1 is x_1, xx_2 is x_2, xx_3 is
% y_1, xx_4 is y_2. The latter two are induced by QSST.
[xx_1,xx_2,xx_3,xx_4] = AutoGillespie_Parallel_Compute(polyfunc1,polyfunc2,Gilltime,initial_condition,NumRealize,1);

% Average over the trajectories
xx1 = sum(xx_1,2)./NumRealize;
xx2 = sum(xx_2,2)./NumRealize;
xx3 = sum(xx_3,2)./NumRealize;
xx4 = sum(xx_4,2)./NumRealize;


figure()
subplot(2,1,2)
hold on
plot(Gilltime,xx1,'b','LineWidth',1.5,'DisplayName','X')
plot(Gilltime,xx2,'r','LineWidth',1.5,'DisplayName','Y')
ylabel('Stochastic','Fontsize',16)
xlabel('Time','Fontsize',16)
legend('X','Y')
subplot(2,1,1)
hold on
plot(Gilltime,y_9(:,1),'b','LineWidth',1.5,'DisplayName','X')
plot(Gilltime,y_9(:,2),'r','LineWidth',1.5,'DisplayName','Y')
legend('X','Y')
title('Detected Data','Fontsize',16)
ylabel('Deterministic','Fontsize',16)

figure()
hold on
plot(xx1,xx2,':b','LineWidth',0.5,'DisplayName','Stochastic Path')
title('Chemical Trajectory','Fontsize',16)
xlabel('Specimen X Molecule Count','Fontsize',16)
ylabel('Specimen Y Molecule Count','Fontsize',16)
legend