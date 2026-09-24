
% poisson2.m  -- solve the Poisson problem u_{xx} + u_{yy} = f(x,y)
% on [ax,bx] x [ay,by].  
% 
% The 5-point Laplacian is used at interior grid points.
% This system of equations is then solved using backslash.
% 
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter3  (2007)

function uvec = poisson5grid(ax, bx, mx, ay, by, my, f, bound)

    hx = (bx-ax)/(mx+1);
    hy = (by-ay)/(my+1);
    x = linspace(ax,bx,mx+2);   % grid points x including boundaries
    y = linspace(ay,by,my+2);   % grid points y including boundaries
    
    [X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
    X = X';                     % transpose so that X(i,j),Y(i,j) are
    Y = Y';                     % coordinates of (i,j) point
    
    Iint = 2:mx+1;              % indices of interior points in x
    Jint = 2:my+1;              % indices of interior points in y
    Xint = X(Iint,Jint);       % interior points
    Yint = Y(Iint,Jint);
    
    rhs = f(Xint,Yint);        % evaluate f at interior points for right hand side
                               % rhs is modified below for boundary conditions.
    
    % set boundary conditions around edges of usoln array:
    
    usoln = bound(X, Y);        % use true solution for this test problem
                                % This sets full array, but only boundary values
                                % are used below.  For a problem where utrue
                                % is not known, would have to set each edge of
                                % usoln to the desired Dirichlet boundary values.
    
    
    % adjust the rhs to include boundary terms:
    rhs(:,1) = rhs(:,1) - usoln(Iint,1)/hy^2;
    rhs(:,my) = rhs(:,my) - usoln(Iint,my+2)/hy^2;
    rhs(1,:) = rhs(1,:) - usoln(1,Jint)/hx^2;
    rhs(mx,:) = rhs(mx,:) - usoln(mx+2,Jint)/hx^2;
    
    
    % convert the 2d grid function rhs into a column vector for rhs of system:
    F = reshape(rhs,mx*my,1);
    
    % form matrix A:
    Ix = speye(mx);
    Iy = speye(my);
    ex = ones(mx, 1);
    ey = ones(my, 1);
    T = spdiags([ex/hx^2 -2*(1/hx^2 + 1/hy^2)*ex ex/hx^2],[-1 0 1],mx,mx);
    S = spdiags([ey ey],[-1 1],my,my)/hy^2;
    A = kron(Iy,T) + kron(S,Ix);
    
    
    % Solve the linear system:
    uvec = A\F;
end

