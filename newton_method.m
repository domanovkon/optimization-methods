function newton_method()
    clc();
    global delta; 
    delta = 1e-6;
    debug = true;
    a = 0;
    b = 1;
    f=@(x) (x-0.777).^4;
    %f=@(x) ((4.*x.^3 + 2.*x.^2 - 4.*x + 2).^(2.^0.5)) + asin((1)./(-x.^2 + x + 5)) - 5.0;
    epsilons = [1e-6];
    for E = epsilons
        [x_star, i, xI] = alg(a, b, f, E, debug);
        xRange = a:0.01:b;
        fprintf('Ньютон\n');
        fprintf('E               N               x*               f*               \n');
        fprintf('%13.6f', E);
        fprintf('%13d', length(xI));
        fprintf('%13.10f', x_star);
        fprintf('%13.6f\n', f(x_star));
        x_star-0.777
        fh = figure('Name', 'E = ' + string(E)); 
        fh.WindowState = 'maximized';
        hold on;
        grid on;
        plot(xRange, f(xRange), 'LineWidth', 1);
        plot(xI, f(xI), 'o', 'LineWidth', 1, 'MarkerSize', 8);
        plot(x_star, f(x_star), '--gs','LineWidth', 3, 'MarkerSize', 8);
        legend('f(x)', '(xi; f(xi))', '(x*; f(x*))');
    end
    fprintf('fminbnd\n');
    fprintf('E               N               x*               f*               \n');
    for E = epsilons
        [x, fval, ~, output] = fminbnd(f, a, b, optimset('TolX', E));
        fprintf('%13.6f', E);
        fprintf('%13d', output.iterations);
        fprintf('%13.10f', x);
        fprintf('%13.6f\n', fval);
        x - 0.777
     end
end

%function df = dTargetFunc(f, x, delta)
%    fxdelta1 = f(x+delta);
%    fxdelta2 = f(x-delta);
%    df = (fxdelta1 - fxdelta2) / 2 / delta;
%end

%function d2f = d2TargetFunc(f, x, delta)
%    fxdelta1 = f(x+delta);
%    fxdelta2 = f(x-delta);
%    d2f = (fxdelta1 - 2 * f(x) + fxdelta2) / delta / delta;
%end

function [df,d2f] = dTargetFunc(f, x)
    global delta;
    fxdelta1 = f(x+delta);
    fxdelta2 = f(x-delta);
    df=(fxdelta1-fxdelta2)/(2*delta);
    d2f=(fxdelta1-2*f(x)+fxdelta2)/(delta^2);
end

function [x, i, xI] = alg(a, b, f, eps, debug)
    x=(a+b)/2;
    h=1;
    i=1;
    xI = [x];
    while ~(abs((x-h)-x)<=eps)
        %delta = eps;
        %delta = 10^-4;
        %f1 = dTargetFunc(f,x, delta);
        %f2 = d2TargetFunc(f,x, delta);
        [f1,f2] = dTargetFunc(f,x);
        dx = x;
        h=f1/f2;
        if(dx-h==Inf || isnan(dx-h))
            break
        end
        x=dx-h;
        if(debug)
            fprintf('i = %d\tx = %.15f\tf(x) = %.15f\n', i, x, f(x))
        end
        xI = [xI x];
        i=i+1;
        [df,ddf] = dTargetFunc(f,x);
        %df = dTargetFunc(f,x, delta);
        %ddf = d2TargetFunc(f,x, delta);
        f1 = double(df);
    end
end