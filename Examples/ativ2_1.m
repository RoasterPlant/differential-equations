function G = Green(x, y, h)    
    if y == 0
        G = 1 - x;
    elseif y == 1
        G = x;
    elseif x <= y
        G = h * (y - 1) * x; 
    else
        G = h * (x - 1) * y;
    end
end

n = 5;
h = 1/(n - 1);
G = @(x, y) Green(x, y, h); % Define Green function for h = 0.25
x = linspace(0, 1, 5);
[X, Y] = meshgrid(x, x);
B = arrayfun(G, X, Y);

plot(x, B(1, :))
xlabel('x');
ylabel('y(x)');
labels = strings(n, 1);
labels(1) = 'G0(x)';
labels(n) = 'G1(x)';
hold on

for i=2:n-1
    plot(x, B(i, :))
    labels(i) = sprintf('G(x; %0.2f)', x(i));
end

plot(x, B(n, :))
legend(labels, 'Location', 'best');

hold off


