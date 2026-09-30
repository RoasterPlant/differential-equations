function f = f(t, u)
    f = [u(2); u(3); 4 * t^2 - u(3) - 4 * u(2) - 4 * u(1)];
end

function u = u_true(t)
    u = [-sin(2 * t) + t.^2 - 3; -2 * (cos(2 * t) - t); 4 * sin(2 * t) + 2];
end

n = 6;
epsilon = logspace(-5, 0, n);
errors = zeros(3, n);
evaluations = zeros(n, 1);
u0 = [-3; -2; 2];
tspan = [0 2];

T = linspace(tspan(1), tspan(2), 1000);
U = u_true(T);

for i = 1:n
    sol = ode113(@f, tspan, u0, odeset('AbsTol', epsilon(i)));
    t = sol.x;
    u = sol.y;
    errors(:, i) = vecnorm(u - u_true(t), inf, 2);
    stats = sol.stats;
    evaluations(i) = stats.nfevals;
end

plot(t, u,'--', T, U, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Solution components');
legend('u_1', 'u_2', 'u_3', 'u_1 analítica', 'u_2 analítica', 'u_3 analítica', 'Location', 'best');
grid on;

% Create error table

errorTable = table(epsilon(:), errors(1, :).', errors(2, :).', errors(3, :).', evaluations, ...
    'VariableNames', {'AbsTol', 'ErrorU1', 'ErrorU2', 'ErrorU3', 'FunctionEvaluations'});
disp(errorTable);