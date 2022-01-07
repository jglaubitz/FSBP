%% script_linear
%
% Description: 
%  Script to numerically solve the linear advection equation 
%  Periodic boundary conditions 
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: Jan 07, 2022

%% Setting up the script 
clc, clear 
 
%% Parameters of the problem 
x_L = 0; x_R = 1; % domain boundaries 
T = 1.0; % end time  
u_init = @(x) cos(4*pi*x) + 0.5*sin(40*pi*x); % initial data 

%% Shared parameters for the SBP-SAT method 
K = 81; % dimension of approximation space 
I = 1; % number of blocks  
x_eval = linspace(x_L,x_R,1000); % evaluation points for reference solution

%% Solve the problem using a polynomial function space on Lobatto points 
approx_space = 'poly'; % approximation space (poly, trig, exp, cubic)  
points = 'Lobatto'; % data points (equid, Lobatto, Halton, random) 
% Solve problem 
[ x_poly, u_poly, mass_poly, energy_poly, u_ref ] = solve_linear_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval );

%% Solve the problem using a trigonometric function space on equidistant points 
approx_space = 'trig'; % approximation space (poly, trig, exp, cubic)  
points = 'equid'; % data points (equid, Lobatto, Halton, random) 
% Solve problem 
[ x_trig, u_trig, mass_trig, energy_trig, u_ref ] = solve_linear_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval );

%% Plots 

% Plot the solutions 
figure(1) 
p = plot( x_eval, u_ref,'k:', x_poly, u_poly,'b--', x_trig, u_trig,'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
ylim([-1.75,1.75]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','trig');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')

% Plot the mass 
t = mass_poly(:,1); % list of times 
%syms x; aux = int(u_init(x),x_L,x_R)
aux = 0;
mass_ref = aux*t.^0; % mass of exact solution over time 
figure(2) 
p = plot( t, mass_ref,'k:', t, mass_poly(:,2),'b--', mass_trig(:,1), mass_trig(:,2),'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$\int u \mathrm{d}x$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','trig');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')

% Plot the energy 
t = energy_poly(:,1); % list of times 
%syms x; aux = int(u_init(x)^2,x_L,x_R)
aux = 5/8;
energy_ref = aux*t.^0; % energy of exact solution over time 
figure(3) 
p = plot( t, energy_ref,'k:', t, energy_poly(:,2),'b--', energy_trig(:,1), energy_trig(:,2),'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$\int u^2 \mathrm{d}x$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','trig');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')