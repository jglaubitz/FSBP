%% characteristic_tracing_Burgers
%
% Description: 
%  Function to compute (smooth) solutions of Burgers' equation
%
% Author: Jan Glaubitz 
% Date: Jan 07, 2022 
% 
% INPUT: 
%  x :          grid points 
%  t :          time
%  u_init :    	initial condition 
%
% OUTPUT: 
%  u :              (approximate) solution values at the grid points
%  

function [u] = characteristic_tracing_Burgers( x, t, u_init )

    [N,I] = size(x); 
    u = zeros(N,I); 

    % root finder - solve the characteristic equation
    for i=1:I
        for n=1:N 
            x0 = u_init(x(n,i));
            fun = @(u) u_init( x(n,i) - t*u ) - u;
            u(n,i) = fzero( fun, x0 ); 
        end
    end

end

