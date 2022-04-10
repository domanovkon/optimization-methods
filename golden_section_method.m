function golden_section_method()
    clc();
    debug = true;
    a = 0;
    b = 1;
    func=@(x) (x-0.777).^14;
    %func=@(x) ((4.*x.^3 + 2.*x.^2 - 4.*x + 2).^(2.^0.5)) + asin((1)./(-x.^2 + x + 5)) - 5.0;
    epsilons = [1e-6];
    for E = epsilons
        [x_star, f_star, xI] = alg(a, b, func, E, debug);
        xRange = a:0.01:b;
        fprintf('E               N               x*               f*               \n');
        fprintf('%13.6f', E);
        fprintf('%13d', length(xI));
        fprintf('%13.10f', x_star);
        fprintf('%13.6f\n', f_star);
        x_star-0.777
        
        fh = figure('Name', 'E = ' + string(E)); 
        fh.WindowState = 'maximized';
        hold on;
        grid on;
        plot(xRange, func(xRange), 'LineWidth', 1);
        plot(xI, func(xI), 'o', 'LineWidth', 1, 'MarkerSize', 8);
        plot(x_star, f_star, '--gs','LineWidth', 3, 'MarkerSize', 8);
        legend('f(x)', '(xi; f(xi))', '(x*; f(x*))');
    end
end

function [x_star, f_star, xI] = alg(a, b, f, eps, debug)
    tau = (sqrt(5) - 1) / 2;
    l = b - a;
    x1 = b - tau * l;
    x2 = a + tau * l;
    f1 = f(x1);
    f2 = f(x2);
    xI = [x1 x2];
    i = 0;
    while l > 2 * eps
        i = i+1;
        if f1 <= f2
            b = x2;
            l = b - a;
            x2 = x1;
            f2 = f1;
            x1 = b - tau * l;
            f1 = f(x1);
            xI = [xI x1];
        else
            a = x1;
            l = b - a;
            x1 = x2;
            f1 = f2;
            x2 = a + tau * l;
            f2 = f(x2);
            xI = [xI x2];
        end
        if (debug)
            fprintf('i = %d\t\tx1 = %.15f\t\tx2 = %.15f\t\tl = %.15f\t|\tf(x1) = %.15f\n', i, x1, x2, l, f(x1))
        end
    end
    x_star = (a + b) / 2;
    f_star = f(x_star);
    xI = [xI x_star];
end
