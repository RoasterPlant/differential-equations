ax = 0;
bx = 1;
ay = 0;
by = 1;
f = @(x,y) 1.25*exp(x+y/2);         % f(x,y) function
bound = @(x,y) exp(x+y/2);

% Solve the problem for nxtest x nytest different grid sizes to test convergence:
mxvals = [10 20 40 80];
myvals = [10 20 40 80];
nxtest = length(mxvals);
nytest = length(myvals);
[mxgrid, mygrid] = meshgrid(mxvals, myvals);
hxvals = (bx - ax)./(mxvals + 1);  % to hold hx values
hyvals = (by - ay)./(myvals + 1); % to hold hy values
E = zeros(nxtest, nytest);      % to hold errors

for i=1:nxtest
    for j=1:nytest
        mx = mxvals(i);
        my = myvals(j);
        uvec = poisson5grid(ax, bx, mx, ay, by, my, f, bound);
        
        x = linspace(ax,bx,mx+2);   % grid points x including boundaries
        y = linspace(ay,by,my+2);   % grid points y including boundaries
        
        [X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
        X = X';                     % transpose so that X(i,j),Y(i,j) are
        Y = Y';                     % coordinates of (i,j) point
        
        utrue = exp(X+Y/2);        % true solution for test problem
        
        % set boundary conditions around edges of usoln array:
        
        usoln = utrue;              % use true solution for this test problem
        % This sets full array, but only boundary values
        % are used below.  For a problem where utrue
        % is not known, would have to set each edge of
        % usoln to the desired Dirichlet boundary values.
        
        % reshape vector solution uvec as a grid function and 
        % insert this interior solution into usoln for plotting purposes:
        % (recall boundary conditions in usoln are already set)
        
        Iint = 2:mx+1;              % indices of interior points in x
        Jint = 2:my+1;              % indices of interior points in y
        %Xint = X(Iint,Jint);       % interior points
        %Yint = Y(Iint,Jint);
        
        usoln(Iint,Jint) = reshape(uvec,mx,my);
        
        % assuming true solution is known and stored in utrue:
        E(i, j) = max(max(abs(usoln-utrue)));
    end
end

% create error table.

fprintf('      hx/hy   ')
fprintf('    %.4e', hyvals)
fprintf('\n\n')
p = [hxvals' E];
for i=1:nxtest
    fprintf('    %.4e', p(i, :))
    fprintf('\n')
end

% fprintf('Error relative to true solution of PDE = %10.3e \n',err)
% 
% % plot results:
% 
% clf
% hold on

% plot grid:
% plot(X,Y,'g');  plot(X',Y','g')

% plot solution:
% contour(X,Y,usoln,50)
% 
% axis([ax bx ay by])
% daspect([1 1 1])
% title('Contour plot of computed solution')
% hold off