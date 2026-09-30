g = @(t) (t * (t * (t * (t - 6) + 12) - 14) + 9) / (1 + t)^2;
f = @(t, u) u^2 - g(t);
t0 = 0;
tf = 1.6;
u0 = 2;
nvals = 2.^(linspace(3, 6, 4));
errors = zeros(size(nvals));
hvals = (tf - t0) ./ nvals;

method_name = "am2"; % Must be one of the following options
switch method_name
    case "eulerIm"
        method = @eulerIm;
    case "eulerEx"
        method = @eulerEx;
    case "rk4"
        method = @rk4;
    case "am2"
        method = @am2;
    otherwise
        error('Unknown method: %s', method_name);
end

u_true = @(t) (1 - t) .* (2 - t) ./ (1 + t);
T = linspace(t0, tf, 1000);

plot(T, u_true(T), 'LineWidth', 1.5)
xlabel('t');
ylabel('u');

hold on

labels = strings(length(nvals) + 1, 1);
labels(1) = 'Analytical solution';

for i=1:length(nvals)
    [t, u] = method(f, [t0 tf], u0, nvals(i));
    scatter(t, u, 10, "filled")
    errors(i) = vecnorm(u - u_true(t), inf, 2);
    labels(i + 1) = sprintf('h = %0.3e', hvals(i));
end

legend(labels, 'Location', 'best');

hold off

error_table(hvals, errors)