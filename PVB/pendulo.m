function T = pendulo(t0, tf, alpha, beta, n, initialGuess, MAX_ITER, epsilon)
    m = n + 1;
    n0 = n - 1;

    % set grid points:  
    gridchoice = 'uniform';          % see xgrid.m for other choices
    t = xgrid(t0, tf, n0, gridchoice);

    switch initialGuess
        case "Linear"
            theta = linearAprox(t, alpha, beta);
        otherwise % case "Zero"
            I = eye(m);
            theta = alpha * I(:, 1) + beta * I(:, m);         
    end
    noBoundTheta = theta(2:end - 1);

    coeff = fdcoeffF(2, t(2), t(1:3));

    for k = 1:MAX_ITER
        J = Jacobian(coeff, noBoundTheta);
        F = - G(coeff, noBoundTheta, alpha, beta);
        d = F.'/J;
        noBoundTheta = noBoundTheta + d.';  % Update the noBoundTheta with the correction
        if max(abs(d)) < epsilon
            break;
        end
    end
    theta(2:end - 1) = noBoundTheta;
    T = theta;  % Store the final approximation in T
end

function A = Jacobian(coeff, theta)
    m = length(theta);
    B = spalloc(m,m,3*m);

    B(1,1:2) = coeff(end - 1:end) + [cos(theta(1)), 0];
    for i = 2:m-1
        B(i, i-1:i+1) = coeff + [0, cos(theta(i)), 0];
    end
    B(m, m-1:m) = coeff(1: 2) + [0, cos(theta(m))];
    A = B;
end

function g = g(coeff, theta, i)
    g = dot(coeff, theta(i-1:i+1)) + sin(theta(i));
end

function G = G(coeff, theta, alpha, beta)
    m = length(theta);
    H = zeros(m, 1);

    H(1) = g(coeff, [alpha, theta(1), theta(2)], 2);
    for i = 2:m-1
        H(i) = g(coeff, theta, i);
    end
    H(m) = g(coeff, [theta(m-1), theta(m), beta], 2);
    G = H;
end

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
        M = [sin(t(1)) cos(t(1))
             sin(t(end)) cos(t(end))];
        c = boundary/M;
    end
    f = @(x) [sin(x), cos(x)] * c.';
    theta = f(t);
end
