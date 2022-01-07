%% script_linear_source
%
% Description: 
%  Script to numerically solve the linear advection equation with a source term. 
%  Inflow boundary conditions. 
%  The FSBP-SAT method is used on a multi-block structure. 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method. 
%
% Author: Jan Glaubitz 
% Date: Jan 07, 2022

%% Setting up the script 
clc, clear 
 
%% Parameters of the problem 
x_L = 0; x_R = pi; % domain boundaries 
T = 3.5; % end time  
source = '2u'; % source term (2u, 2xu)
u_init = @(x) x.^0; % initial data 

%% Shared parameters for the SBP-SAT method 
K = 3; % dimension of approximation space 
I = 6; % number of blocks  
x_eval = linspace(x_L,x_R,1000); % evaluation points for reference solution

%% Solve the problem using a polynomial function space on Lobatto points 
approx_space = 'poly'; % approximation space (poly, trig, exp, cubic)  
points = 'Lobatto'; % data points (equid, Lobatto, Halton, random) 
% Solve problem 
[ x_poly, u_poly, u_ref ] = solve_linear_source_SAT( x_L, x_R, T, source, u_init, I, approx_space, K, points, x_eval );

%% Solve the problem using a trigonometric function space on equidistant points 
approx_space = 'exp'; % approximation space (poly, trig, exp, cubic)  
points = 'equid'; % data points (equid, Lobatto, Halton, random) 
% Solve problem 
[ x_exp, u_exp, u_ref ] = solve_linear_source_SAT( x_L, x_R, T, source, u_init, I, approx_space, K, points, x_eval );

%% Plot the solutions 
figure(1) 
p = plot( x_eval(:), u_ref(:),'k:', x_poly(:), u_poly(:),'b--', x_exp(:), u_exp(:),'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','exp');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best') 