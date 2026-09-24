t0 = 0;
tf = 27;
alpha = 0.7;  % Dirichlet boundary condition at t0
beta = 0.7;   % Dirichlet boundary condition at tf
MAX_ITER = 50;
epsilon = 1e-5;
initialGuess = "Zero";

% Solve the problem for ntest different grid sizes to test convergence:
m1vals = [10 20 40 80 140 200];
ntest = length(m1vals);
hvals = zeros(ntest,1);  % to hold h values

for jtest=1:ntest
    m1 = m1vals(jtest);
    m2 = m1 + 1;
    m = m1 - 1;                 % number of interior grid points

    % set grid points:  
    gridchoice = 'uniform';          % see xgrid.m for other choices
    t = xgrid(t0,tf,m,gridchoice); 

    hvals(jtest) = (tf-t0)/m1;  % average grid spacing, for convergence tests
    theta = pendulo(t0, tf, alpha, beta, m1, initialGuess, MAX_ITER, epsilon);

    clf
    plot(t,theta,'o')  % plot computed solution
    title(sprintf('Computed solution with %i grid points',m2));

    % pause to see this plot:  
    drawnow
    input('Hit <return> to continue ');

end