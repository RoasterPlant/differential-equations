%
% bvp_2.m 
% second order finite difference method for the bvp
%   u''(x) + u(x) = 0,   u(t0)=alpha,   u(tf)=beta
% Using 3-pt differences on an arbitrary nonuniform grid.
% Should be 2nd order accurate if grid points vary smoothly, but may
% degenerate to "first order" on random or nonsmooth grids. 

t0 = 0;
tf = pi;
alpha = 2;   % Dirichlet boundary condition at t0
beta = 3;    % Dirichlet boundary condition at tf

utrue = @(x) linearAprox(x, alpha, beta);  % true solution

% true solution on fine grid for plotting:
xfine = linspace(t0,tf,101);
ufine = utrue(xfine');

% Solve the problem for ntest different grid sizes to test convergence:
nvals = [10 20 40 80];
ntest = length(nvals);
hvals = zeros(ntest,1);  % to hold h values
E = zeros(ntest,1);      % to hold errors

for jtest=1:ntest
    n = nvals(jtest);
    m2 = n + 1;
    m = n - 1;                 % number of interior grid points

    % set grid points:  
    gridchoice = 'uniform';          % see xgrid.m for other choices
    t = xgrid(t0,tf,m,gridchoice); 

    hvals(jtest) = (tf-t0)/n;  % average grid spacing, for convergence tests
    U = linearPendulum(t0, tf, alpha, beta, n);

    % compute error at grid points:
    uhat = utrue(t);
    err = U - uhat;
    E(jtest) = max(abs(err));  
    disp(' ')
    fprintf('Error with %i points is %9.5e\n',m2,E(jtest))

    clf
    scatter(t,U, 30)  % plot computed solution  % plot computed solution
    xlabel('t')
    ylabel('Theta')
    hold on
    plot(xfine,ufine)  % plot true solution
    legend(sprintf('%i pontos na malha', m2),'Solução Analítica','Location','best')
    hold off

    % pause to see this plot:  
    drawnow
    input('Hit <return> to continue ');

end

error_table(hvals, E);   % print tables of errors and ratios
error_loglog(hvals, E);  % produce log-log plot of errors and least squares fit