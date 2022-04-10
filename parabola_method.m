function parabola_method()
    clc();
    debug = true;
    a = 0;
    b = 1;
    fprintf('E               N               x*               f*               \n');
    f=@(x) (x-0.777).^8;
    %f=@(x) ((4.*x.^3 + 2.*x.^2 - 4.*x + 2).^(2.^0.5)) + asin((1)./(-x.^2 + x + 5)) - 5.0;
    epsilons = [1e-6];
    for E = epsilons
        [x_star, f_star, xI, i] = alg(a, b, f, E, debug);
        xRange = a:0.01:b;
        fprintf('%13.6f', E);
        fprintf('%13d', i);
        fprintf('%13.10f', x_star);
        fprintf('%13.6f\n', f_star);
        x_star-0.777
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

function [xRes, yRes, xCalc, i] = alg(a, b, f, eps, debug)
	[x1, y1, x2, y2, x3, y3] = getBasePoints(a, b, f);
    xCalc = [x1 x2 x3];
    is_first = 1;
    i = 0;
    while 1
        xAvg = getXAvg(x1, y1, x2, y2, x3, y3);
        if is_first
            is_first = 0;
            i = i +1;
        else
            i = i +1;
            if (debug)
                fprintf('i = %d\tx_min = %.15f\tf(x_min) = %.15f\n', i, xAvg, f(xAvg))
            end
            if abs(xAvg - xPrev) <= eps
                xRes = xAvg;
                yRes = f(xAvg);
                xCalc = [xCalc xAvg];
                break;
            end
        end
        xPrev = xAvg;
        yAvg = f(xAvg);
        xCalc = [xCalc xAvg];
        [x1, y1, x2, y2, x3, y3] = getNextPoints(x1, y1, x2, y2, x3, y3, xAvg, yAvg);
    end
end

function [x1, y1, x2, y2, x3, y3] = getBasePoints(a, b, f)
    % Метод золотого сечения
    tau = (sqrt(5) - 1) / 2;
    l = b - a;
    x1 = b - tau * l;
    x2 = a + tau * l;
    f1 = f(x1);
    f2 = f(x2);
    if f2 > f1   
        % Отрезок [a; x2].
        [x1, y1, x2, y2, x3, y3] = deal(a, f(a), x1, f1, x2, f2);
    else
        % Отрезок [x1; b].
        [x1, y1, x2, y2, x3, y3] = deal(x1, f1, x2, f2, b, f(b));
    end
end

function xAvg = getXAvg(x1, f1, x2, f2, x3, f3)
    % Вычисление точки минимума xAvg полинома
	a1 = (f2 - f1) / (x2 - x1);
    a2 = ((f3 - f1) / (x3 - x1) - a1) / (x3 - x2);
	xAvg = (x1 + x2 - a1 / a2) / 2;
end

function [x1, y1, x2, y2, x3, y3] = getNextPoints(x_1, y_1, x_2, y_2, x_3, y_3, xAvg, yAvg)
    % Выбор x1, x2, x3 для остальных итераций с помощью метода исключения отрезков.
    if x_2 > xAvg
        if y_2 > yAvg
            [x1, y1, x2, y2, x3, y3] = deal(x_1, y_1, xAvg, yAvg, x_2, y_2);
        else
            [x1, y1, x2, y2, x3, y3] = deal(xAvg, yAvg, x_2, y_2, x_3, y_3);
        end
    else
        if yAvg > y_2
            [x1, y1, x2, y2, x3, y3] = deal(x_1, y_1, x_2, y_2, xAvg, yAvg);
        else
            [x1, y1, x2, y2, x3, y3] = deal(x_2, y_2, xAvg, yAvg, x_3, y_3);
        end
    end
end