%% solve_Burgers_SAT
%
% Description: 
%  Function to numerically solve the nonlinear Burgers' equation. 
%  Periodic boundary conditions 
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: Jan 07, 2022
% 
% INPUT: 
%  x_L, x_R :       left and right boundary of the domain 
%  T :              end time 
%  I :              number of blocks 
%  approx_space :  	approximation space (poly, trig, exp, cubic) 
%  K :              dimension of the approximation space 
%  points           data points (equid, Lobatto, Halton, random) 
%  x_eval :         points at which the reference solution is evaluated 
%
%
% OUTPUT: 
%  x :      grid points 
%  u_num :  numerical solution at grid points 
%  u_ref :  reference solution    

function [ x, u_num, u_ref ] = solve_Burgers_SAT( x_L, x_R, T, u_init, I, approx_space, K, points, x_eval )

    %% Set-up the method 

    % Data points and the FSBP operator on the reference block [0,1]
    [ x_ref, w_ref ] = compute_QF( 0, 1, approx_space, points, K ); % grid points and weights on the reference block
    N = length(x_ref); % number of data points
    block_width = (x_R-x_L)/I; % block width
    [ basis_F, dx_basis_F, span_G, m_G ] = generate_span( 0, 1, approx_space, points, K ); % bases of different spaces 
    [D, P, Q] = compute_FSBP( basis_F, dx_basis_F, x_ref, w_ref ); % FSBP operator 
    D = (1/block_width)*D; P = block_width*P; 
    P_inv = sparse(inv(P)); % precompute inverse diagonal-norm matrix 
    
    % Time step size
    dx_min = min(x_ref(2:end)-x_ref(1:end-1)); % minimum distance between any two neighboring grid points 
    dt = 0.01*(dx_min*block_width); % time-step size
    
    % Global grid points 
    x = zeros(N,I); 
    for i=1:I 
        x(:,i) = x_L + (x_R-x_L)*(i-1)/I + x_ref*block_width;
    end

    % initial data 
    u = u_init(x); % solution values on the global grid
    
    % reference solution 
    if x_eval==0
        x_eval = x; 
    end
    u_ref = characteristic_tracing_Burgers( x_eval, T, u_init );


    %% Iterate over time with a 3th-order Runge-Kutta until time T is reached 
    t = 0; 
    while (t<T)  

        % time stepping 
        if T-t<dt 
            dt = T-t; 
        else
            t = t+dt
        end
            
      	% 1st update step 
        SAT = compute_SAT( u );
      	for i = 1:I 
           	k1(:,i) = u(:,i) + dt*( -D*u(:,i).^2/3 - u(:,i).*(D*u(:,i))/3 + P_inv*SAT(:,i) );  
        end
        % 2nd update step 
        SAT = compute_SAT( k1 );
        for i = 1:I 
            k1(:,i) = (3/4)*u(:,i) + (1/4)*k1(:,i) + (1/4)*dt*( -D*k1(:,i).^2/3 - k1(:,i).*(D*k1(:,i))/3 + P_inv*SAT(:,i) );  
        end    
        % 3th update step 
        SAT = compute_SAT( k1 ); 
        for i = 1:I
            u(:,i) = (1/3)*u(:,i) + (2/3)*k1(:,i) + (2/3)*dt*( -D*k1(:,i).^2/3 - k1(:,i).*(D*k1(:,i))/3 + P_inv*SAT(:,i) );
        end
        
    end

    u_num = u;
    
end


%% Function to compute BCs and SATs 
function [ SAT ] = compute_SAT( u ) 
    
    [ N, I ] = size(u);
    SAT = zeros(N,I); % initialize SAT
    
    for i=1:I 
        % left boundary condition 
        if i==1 
            g_L = u(N,I); 
        else 
            g_L = u(N,i-1); 
        end
        % right boundary condition 
        if i==I 
            g_R = u(1,1); 
        else 
            g_R = u(1,i+1); 
        end
        % SATs 
        SAT(1,i) = -(2/3)*u(1,i)*( u(1,i) - g_L ); 
        
    end
    
end