function [t, y] = RungeKutta(t0, y0, h, tn, A);
    t = (t0:h:tn)';
    s = size(A);
    y = zeros(s(1), length(t));
    
    y(:,1) = y0;
    for i = 1:1:length(t)-1
        k1 = A * y(:,i);
        k2 = A * (y(:,i) + k1 * h / 2);
        k3 = A * (y(:,i) + k2 * h / 2);
        k4 = A * (y(:,i) + k3 * h);
        y(:,i+1) = y(:,i) + h * (k1 + 2 * k2 + 2 * k3 + k4) /6;
    end
end
