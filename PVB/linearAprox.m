% This is the exact approximation to the linear pendulum problem.
% It MUST be on a separate file because both pendulum.m and
% ativ2_6.m use it, and MATLAB does not allow multiple public functions
% to be defined on a single file.
% It's a stupid feature.

function theta = linearAprox(t, alpha, beta)
    periodicity = (t(end) - t(1))/(2 * pi);
    if isapprox(periodicity,round(periodicity))
        if isapprox(sin(t(1)), 0)
            c = [alpha, alpha];
        else
            c = [alpha/sin(t(1)), 0];
        end
    else
        boundary = [alpha, beta];
        M = [sin(t(1)) sin(t(end)) 
            cos(t(1)) cos(t(end))];
        c = boundary/M;
    end
    f = @(x) [sin(x), cos(x)] * c.';
    theta = f(t);
end