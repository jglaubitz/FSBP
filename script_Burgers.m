%% script_Burgers 
%
% Description: 
%  Script to numerically solve the Burgers equation 
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
T = 0.01; % end time  
u_init = @(x) 1 + 0.5*sin(4*pi*x).^3 + 0.25*cos(4*pi*x).^5; % initial data 

%% Shared parameters for the SBP-SAT method 
K = 3; % dimension of approximation space 
I = 10; % number of blocks  
x_eval = linspace(x_L,x_R,1000); % evaluation points for reference solution

%% Solve the problem using a polynomial function space on Lobatto points 
approx_space = 'poly'; % approximation space (poly, trig, exp, cubic)  
points = 'Lobatto'; % data points (equid, Lobatto, Halton, random) 
% Solve problem 
[ x_poly, u_poly, u_ref ] = solve_Burgers_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval );

%% Solve the problem using a trigonometric function space on equidistant points 
approx_space = 'exp'; % approximation space (poly, trig, exp, cubic)  
points = 'equid'; % data points (equid, Lobatto, Halton, random) 
% Solve problem 
[ x_exp, u_exp, u_ref ] = solve_Burgers_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval );

%% Plot the solutions 
figure(1) 
p = plot( x_eval(:), u_ref(:),'k:', x_poly(:), u_poly(:),'b--', x_exp(:), u_exp(:),'r-' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
ylim([0.4,1.6]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','exp');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')
