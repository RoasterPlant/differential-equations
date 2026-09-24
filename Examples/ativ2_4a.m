%
% bvp_2.m 
% second order finite difference method for the bvp
%   u''(x) = f(x),   u'(ax)=alpha,   u(bx)=sigma
% Using 3-pt differences on an arbitrary nonuniform grid.
% Should be 2nd order accurate if grid points vary smoothly, but may
% degenerate to "first order" on random or nonsmooth grids. 
%
% Different BCs can be specified by changing the first and/or last rows of 
% A and F.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

ax = 0;
bx = 3;
alpha = -5; ao = 0;   % Dirichlet boundary condition at ax
sigma = 3; bo = 1;    % Neumann boundary condition at bx

f = @(x) exp(x);  % right hand side function
utrue = @(x) exp(x) + (sigma - exp(bx))*(x - ax) + alpha - exp(ax);  % true solution

% true solution on fine grid for plotting:
xfine = linspace(ax,bx,101);
ufine = utrue(xfine);

% Solve the problem for ntest different grid sizes to test convergence:
m1vals = [10 20 40 80];
ntest = length(m1vals);
hvals = zeros(ntest,1);  % to hold h values
E = zeros(ntest,1);      % to hold errors

for jtest=1:ntest
    m1 = m1vals(jtest);
    m2 = m1 + 1;
    m = m1 - 1;                 % number of interior grid points

    % set grid points:  
    gridchoice = 'rtlayer';          % see xgrid.m for other choices
    x = xgrid(ax,bx,m,gridchoice); 

    hvals(jtest) = (bx-ax)/m1;  % average grid spacing, for convergence tests
    U = bvp2(ax, bx, alpha, ao, sigma, bo, f, m1);

    % compute error at grid points:
    uhat = utrue(x);
    err = U - uhat;
    E(jtest) = max(abs(err));  
    disp(' ')
    fprintf('Error with %i points is %9.5e\n',m2,E(jtest))

    clf
    plot(x,U,'o')  % plot computed solution
    title(sprintf('Computed solution with %i grid points',m2));
    hold on
    plot(xfine,ufine)  % plot true solution
    hold off

    % pause to see this plot:  
    drawnow
    input('Hit <return> to continue ');

end

error_table(hvals, E);   % print tables of errors and ratios
error_loglog(hvals, E);  % produce log-log plot of errors and least squares fit