function [t, u] = am2(f, tspan, u0, steps)
    t = linspace(tspan(1), tspan(2), steps + 1);
    u = zeros(numel(u0), steps + 1);
    fu = zeros(numel(u0), steps + 1);
    u(:, 1) = u0;
    fu(:, 1) = f(t(1), u0);
    stepSize = (tspan(2) - tspan(1)) / steps;

    % Initial step (trapezoidal method)
    trap = @(next) next - u(:, 1) - stepSize/2 * (f(t(2), next) + fu(1));
    u(:, 2) = fsolve(trap, u(:, 1), optimoptions("fsolve", Display="none"));
    fu(:, 2) = f(t(2), u(:, 2));

    % 2-step Adam-Moulton
    aux = @(u_next, u_prev, fu, t) (u_next - u_prev) - stepSize/12 * (-fu(1) + 8 * fu(2) + 5 * f(t, u_next)); 
    for i = 1:steps - 1
        F = @(u_next) aux(u_next, u(:, i + 1), fu(:, i:i + 1), t(i + 2));
        u(:, i + 2) = fsolve(F, u(:, i + 1), optimoptions("fsolve", Display="none"));
        fu(:, i + 2) = f(t(i + 2), u(:, i + 2));
    end
end