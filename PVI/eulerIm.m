function [t, u] = eulerIm(f, tspan, u0, steps)
    t = linspace(tspan(1), tspan(2), steps + 1);
    u = zeros(numel(u0), steps + 1);
    u(:, 1) = u0;
    stepSize = (tspan(2) - tspan(1)) / steps;
    aux = @(u_next, u_prev, t) (u_next - u_prev) - stepSize * f(t, u_next); 
    for i = 1:steps
        F = @(u_next) aux(u_next, u(:, i), t(i));
        u(:, i + 1) = fsolve(F, u(:, i), optimoptions("fsolve", Display="none"));
    end
end