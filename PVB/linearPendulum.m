% Finite difference solution to the linear problem
function T = linearPendulum(t0, tf, alpha, beta, n)
    m = n + 1;
    n0 = n - 1;                 % number of interior grid points
    
    % set grid points:  
    gridchoice = 'uniform';          % see xgrid.m for other choices
    x = xgrid(t0,tf,n0,gridchoice);   
    
    % set up matrix A (using sparse matrix storage):
    A = spalloc(m,m,3*m);   % initialize to zero matrix
    
    % first row:
    A(1,1:3) = fdcoeffF(0, x(1), x(1:3)); 
    
    % interior rows:
    for i=2:n
        A(i,i-1:i+1) = fdcoeffF(2, x(i), x((i-1):(i+1)));
    end
    
    % last row:
    A(m,n0:m) = fdcoeffF(0, x(m), x(n0:m));
    
    % Add u(x):
    A(2:n, 2:n) = A(2:n, 2:n) + eye(n0);
    
    % Right hand side:
    F = zeros(m, 1); 
    F(1) = alpha;  
    F(m) = beta;
    
    % solve linear system:
    T = A\F;
end