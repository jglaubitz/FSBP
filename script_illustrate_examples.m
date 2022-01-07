%% script_illustrate_examples
%
% Description: 
% Script to illustrate the example's solutions 
%
% Author: Jan Glaubitz 
% Date: Jan 07, 2022

clear all; close all; clc; % clean up
%warning('off','all') % in case any of the warnings become too anoying 


%% Set up grid points 
n = 1000; % number of grid points 
x = linspace(0,1,n)'; % equidistant grid points 


%% Boundary layer problem 
eps = 10^(-2); norm_const = exp(1/(2*eps)) - 1;
fun_BL = @(x) ( exp(x/(2*eps)) - 1 )/norm_const; % Function  
y_BL = fun_BL(x); % function values at the grid points 


%% Highly oscillatory function  
fun_osc = @(x) cos(4*pi*x) + 0.5*sin(40*pi*x);  % Function
y_osc = fun_osc(x); % function values at the grid points 
plot(x,y_osc)


%% Illustrate for boundary layer function (k=3)

% Polynomial of degree at most 2 
[ fit_poly_BL, gof, fitinfo ] = fit( x, y_BL, 'poly2'); 
%res = fitinfo.residuals; % vector of residuals 
%L2_error_poly_osc = sqrt( norm(res)/norm(y_osc) ); % L^2-error 
%Linf_error_poly_osc = max(res)/max(y_osc); % L^inf-error 

% Exponential space 
ft = fittype('a + b*x + c*exp(0.5*(10^2)*x)');
[ fit_exp_BL, gof, fitinfo ] = fit( x, y_BL, ft, 'StartPoint',[-1/norm_const,0,1/norm_const] ); 

% Plot the results 
figure(1) 
p1 = fplot( fun_BL, [0,1], 'k:' ); % plot the reference solution 
set(p1, 'LineWidth',3); 
hold on 
p2 = plot( x, fit_poly_BL(x), 'b--' ); 
set(p2, 'LineWidth',3); 
p3 = plot( x, fit_exp_BL(x), 'r-.' ); 
set(p3, 'LineWidth',3); 
set(gca, 'FontSize', 24); % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex'); 
ylabel('$u$','Interpreter','latex'); 
lgnd = legend('ref','poly','exp','Interpreter','latex','Location','best');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on 
hold off 


%% Illustrate for highly oscillatory function (k=3)

% Polynomial of degree at most 2 
[ fit_poly_osc, gof, fitinfo ] = fit( x, y_osc, 'poly2'); 

% Trigonometric space ( alpha=4*pi, xi=0 ) 
ft = fittype('a + b*sin( 4*pi*x ) + c*cos( 4*pi*x )');
[ fit_tri_osc, gof, fitinfo ] = fit( x, y_osc, ft, 'StartPoint',[0,0,1] ); 

% Plot the results 
figure(2) 
p1 = fplot( fun_osc, [0,1], 'k:' ); % plot the reference solution 
set(p1, 'LineWidth',3); 
hold on 
p2 = plot( x, fit_poly_osc(x), 'b--' ); 
set(p2, 'LineWidth',3); 
p3 = plot( x, fit_tri_osc(x), 'r-.' ); 
set(p3, 'LineWidth',3); 
set(gca, 'FontSize', 24); % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex'); 
ylabel('$u$','Interpreter','latex'); 
lgnd = legend('ref','poly','trig','Interpreter','latex','Location','best');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on 
hold off 