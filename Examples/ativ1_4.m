
% Example 1.1 
% Table 1.1 and Figure 1.2 of  http://www.amath.washington.edu/~rjl/fdmbook

function d = Dplus(u, x, h)
    d = (u(x+h) - u(1))/h;
end

function d = Dzero(u, x, h)
d = (u(x + h) - u(x - h))/(2*h);
end

function d = D3(u, x, h)
d = (2*u(x+h) + 3*u(x) - 6*u(x - h) + u(x - 2*h))/(6*h);
end

hvals = logspace(-4, -1, 20);
Epu = []; E0u = []; E3u = [];
% u = @(x) sin(10*x);
% uptrue = 10 * cos(10);
% u = @(x) x*(x*(4*x - 12) + 14) - 4;
% uptrue = x*(12*x - 24) + 14;
u = @(x) abs(x - 1.001).^(3/2);
uptrue = - (3/2)*0.001.^(1/2);
x0 = 1;

% table headings:
disp(' ')
disp('       h              Dpu             D0u             D3u')

for i=1:length(hvals)
    h = hvals(i);
    % approximations to u'(1):
    Dpu = Dplus(u, x0, h);
    D0u = Dzero(u, x0, h);
    D3u = D3(u, x0, h);
    % errors:
    Epu(i) = Dpu - uptrue;
    E0u(i) = D0u - uptrue;
    E3u(i) = D3u - uptrue;

    % print line of table:
    fprintf('%13.4e   %13.4e   %13.4e  %13.4e\n',...
        h,Epu(i),E0u(i),E3u(i))
end

% plot errors:
clf
loglog(hvals,abs(Epu),'o-')
% axis([5e-4 .2 1e-12 1])
hold on
loglog(hvals,abs(E0u),'o-')
loglog(hvals,abs(E3u),'o-')
hold off

