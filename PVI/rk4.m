function [t, u] = rk4(f, tspan, u0, steps)
    t = linspace(tspan(1), tspan(2), steps + 1);
    u = zeros(numel(u0), steps + 1);
    u(:, 1) = u0;
    stepSize = (tspan(2) - tspan(1)) / steps;
    for i = 1:steps
        k1 = f(t(i), u(:, i));
        k2 = f(t(i) + stepSize / 2, u(:, i) + stepSize * k1 / 2);
        k3 = f(t(i) + stepSize / 2, u(:, i) + stepSize * k2 / 2);
        k4 = f(t(i) + stepSize, u(:, i) + stepSize * k3);
        u(:, i + 1) = u(:, i) + stepSize * (k1 + 2 * k2 + 2 * k3 + k4)/6;
    end
end