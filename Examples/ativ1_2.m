function y = f(x)
    y = sin(2 * x);
end

function y = df(x)
    y =-4 * sin(2 * x);
end

hvals = logspace(-1, -4, 13);
[coefficients, error] = fdstencil(2, -2:2);
uptrue = df(1);

E2f = [];

% table headings:
disp(' ')
disp('       h              Error        Pred_Error')

for i=1:length(hvals)
    h = hvals(i);
    % approximations to f'(1):
    D2f = 0;
    for j=-2:2
        D2f = D2f + coefficients(j + 3) * f(1 + j * h) / h.^2;
    end
    % errors:
    E2f(i) = abs(D2f - uptrue);
    E2f_pred(i) = abs(error(1) * h.^(error(2)) * 32 * cos(2) + error(4) * h.^(error(5)) * (-64) * f(1));

    % print line of table:
    fprintf('%13.4e   %13.4e   %13.4e\n',...
        h,E2f(i), E2f_pred(i))
end

% plot errors:
clf
loglog(hvals, E2f,'o-')
axis([1e-4 .1 1e-12 1])