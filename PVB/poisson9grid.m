
% poisson2.m  -- solve the Poisson problem u_{xx} + u_{yy} = f(x,y)
% on [ax,bx] x [ay,by].  
% 
% The 5-point Laplacian is used at interior grid points.
% This system of equations is then solved using backslash.
% 
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter3  (2007)

function uvec = poisson9grid(ax, bx, mx, ay, by, my, f, bound)
    
    hx = (bx-ax)/(mx+1);
    hy = (by-ay)/(my+1);
    x = linspace(ax,bx,mx+2);   % grid points x including boundaries
    y = linspace(ay,by,my+2);   % grid points y including boundaries
    
    [X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
    X = X';                     % transpose so that X(i,j),Y(i,j) are
    Y = Y';                     % coordinates of (i,j) point

    rhs = f(X,Y);        % evaluate f at interior points for right hand side
    rhs = stencil5(rhs);        % Add 2-order error
    
    A = 12 * (hx * hy)^2 / (hx^2 + hy^2); % 6 h^2
    B = (10 * hy^2 - 2 * hx^2) / (hx^2 + hy^2); % 4
    C = (10 * hx^2 - 2 * hy^2) / (hx^2 + hy^2); % 4
    
    % rhs is modified below for boundary conditions.
    
    % set boundary conditions around edges of usoln array:
    
    usoln = bound(X, Y);
    % This sets full array, but only boundary values
    % are used below.  For a problem where utrue
    % is not known, would have to set each edge of
    % usoln to the desired Dirichlet boundary values.
    
    
    % adjust the rhs to include boundary terms:
    rhs(:,1) = rhs(:,1) - convn(usoln(:, 1), [1; C; 1], "valid")/A;
    rhs(:,my) = rhs(:,my) - convn(usoln(:, my + 2), [1; C; 1], "valid")/A;
    rhs(1,:) = rhs(1,:) - convn(usoln(1, 2:my+1), [1 B 1], "same")/A;
    rhs(mx,:) = rhs(mx,:) - convn(usoln(mx + 2, 2:my+1), [1 B 1], "same")/A;
    
    
    % convert the 2d grid function rhs into a column vector for rhs of system:
    F = reshape(rhs,mx*my,1);
    
    % form matrix G for discrete linear system:
    I = speye(my);
    ex = ones(mx, 1);
    ey = ones(my, 1);
    T = spdiags([B * ex -20 * ex B * ex],[-1 0 1],mx,mx);
    S = spdiags([ey ey],[-1 1],my,my);
    R = spdiags([ex C * ex ex],[-1 0 1],mx,mx);
    G = (kron(I,T) + kron(S,R)) / A;
    
    
    % Solve the linear system:
    uvec = G\F;
end

function st = stencil5(x)
    A = 1/12 *[0 1 0;
               1 8 1;
               0 1 0];
    st = conv2(x, A, "valid");
end
