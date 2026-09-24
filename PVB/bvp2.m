%
% bvp_2.m 
% second order finite difference method for the bvp
%   u''(x) = f(x),   u^(ao)(ax)=alpha,   u^(bo)(bx)=beta
% Using 3-pt differences on an arbitrary nonuniform grid.
% Should be 2nd order accurate if grid points vary smoothly, but may
% degenerate to "first order" on random or nonsmooth grids. 
%
% Different BCs can be specified by changing the first and/or last rows of 
% A and F.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

function U = bvp2(ax, bx, alpha, ao, beta, bo, f, m1)
    m2 = m1 + 1;
    m = m1 - 1;                 % number of interior grid points

    % set grid points:  
    gridchoice = 'rtlayer';          % see xgrid.m for other choices
    x = xgrid(ax,bx,m,gridchoice);   

    % set up matrix A (using sparse matrix storage):
    A = spalloc(m2,m2,3*m2);   % initialize to zero matrix

    % first row:
    A(1,1:3) = fdcoeffF(ao, x(1), x(1:3)); 

    % interior rows:
    for i=2:m1
        A(i,i-1:i+1) = fdcoeffF(2, x(i), x((i-1):(i+1)));
    end

    % last row:
    A(m2,m:m2) = fdcoeffF(bo, x(m2), x(m:m2)); 

    % Right hand side:
    F = f(x); 
    F(1) = alpha;  
    F(m2) = beta;

    % solve linear system:
    U = A\F;
end
