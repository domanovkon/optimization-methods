function bitwise_search_method123()
    clc();
    debug = true;
    %f=@(x) ((4.*x.^3 + 2.*x.^2 - 4.*x + 2).^(2.^0.5)) + asin((1)./(-x.^2 + x + 5)) - 5.0;
    f=@(x) (x-0.777).^4;
    a = 0;
    b = 1;
    epsilons = [1e-4];
    for E = epsilons
        [x_star, f_star, xI] = alg(a, b, f, E, debug);
        xRange = a:0.01:b;
        fh = figure('Name', 'E = ' + string(E)); 
        fh.WindowState = 'maximized';
        hold on;
        grid on;
        plot(xRange, f(xRange), 'LineWidth', 1);
        plot(xI, f(xI), 'o', 'LineWidth', 1, 'MarkerSize', 8);
        plot(x_star, f_star, '--gs','LineWidth', 3, 'MarkerSize', 8);
        legend('f(x)', '(xi; f(xi))', '(x*; f(x*))');
    end
end

function [x_star, f_star, xI] = alg(a, b, f, eps, debug)
    delta=(b-a)/4;
    x0 = a;
    i = -1;
    xI = [x0];
    
    while true
        x1 = x0 + delta;
        i = i+1;
        
        if(debug)
            fprintf('i = %d\t\tx0 = %.10f\t\tx1 = %.10f\t\tdelta = %.10f\t\tf(x0) = %.10f\n', i, x0, x1, delta, f(x0))
        end
        xI = [xI x1];
        if f(x0)>f(x1)
            x0=x1;
            
            if ~((a<x0) && (x0<b))
                if abs(delta)<=eps
                    break
                else
                    x0=x1;
                    delta=-delta/4;
                end
            end
        elseif abs(delta)<=eps
                    break
        else
            x0=x1;
            delta=-delta/4;
        end
    end
    x_star = x0;
    f_star = f(x_star);
    x_star - 0.777
    fprintf('Итераций = %d\t\txmin = %.10f\t\tfmin = %.10f\n\n', i, x0, f(x0))
end
