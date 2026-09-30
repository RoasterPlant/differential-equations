function [t, u] = eulerEx(f, tspan, u0, steps)
    t = linspace(tspan(1), tspan(2), steps + 1);
    u = zeros(numel(u0), steps + 1);
    u(:, 1) = u0;
    stepSize = (tspan(2) - tspan(1)) / steps;
    for i = 1:steps
        u(:, i + 1) = u(:, i) + stepSize * f(t(i), u(:, i));
    end
end