%
% bvp4.m 
% second order finite difference method for the bvp
%   u''(x) = f(x),   u'(ax)=alpha,   u(bx)=beta
% fourth order finite difference method for the bvp
%   u'' = f,   u'(ax)=alpha,   u(bx)=beta
% Using 5-pt differences on an arbitrary grid.
% Should be 4th order accurate if grid points vary smoothly.
%
% Different BCs can be specified by changing the first and/or last rows of 
% A and F.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter2  (2007)

function U = bvp4(ax, bx, alpha, ao, beta, bo, f, m1)
  m2 = m1 + 1;
  m = m1 - 1;                 % number of interior grid points

  % set grid points:  
  gridchoice = 'uniform';
  x = xgrid(ax,bx,m,gridchoice);   

  % set up matrix A (using sparse matrix storage):
  A = spalloc(m2,m2,5*m2);   % initialize to zero matrix

  % first row for Neumann BC on u'(x(1))
  A(1,1:5) = fdcoeffF(ao, x(1), x(1:5));
  % second row for u''(x(2))
  A(2,1:6) = fdcoeffF(2, x(2), x(1:6));

  % interior rows:
  for i=3:m
     A(i,i-2:i+2) = fdcoeffF(2, x(i), x((i-2):(i+2)));
  end

  % next to last row for u''(x(m+1))
  A(m1,m-3:m2) = fdcoeffF(2,x(m1),x(m-3:m2));
  % last row for Dirichlet BC on u(x(m+2))
  A(m2,m-2:m2) = fdcoeffF(bo,x(m2),x(m-2:m2));

  % Right hand side:
  F = f(x); 
  F(1) = alpha;  
  F(m2) = beta;
  
  % solve linear system:
  U = A\F;

end