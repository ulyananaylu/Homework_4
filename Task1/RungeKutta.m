function [t, y] = RungeKutta(t0, y0, h, tn, A)
    fprintf('Runge-Kutta\n\n');
    fprintf('%10s%30s%30s%30s%30s%30s\n', 'i', 'y1i', 'y2i', 'ti', 'dy1/dt(yi,ti)', 'dy2/dt(yi,ti)');
    dy = A * y0;
    fprintf('%10d%30.4f%30.4f%30.4f%30.4f%30.4f\n', 0, y0(1), y0(2), t0, dy(1), dy(2));
    t = (t0:h:tn)';
    y = zeros(2, length(t));
    
    y(:,1) = y0;
    for i = 1:1:length(t)-1
        k1 = A * y(:,i);
        k2 = A * (y(:,i) + k1 * h / 2);
        k3 = A * (y(:,i) + k2 * h / 2);
        k4 = A * (y(:,i) + k3 * h);
        y(:,i+1) = y(:,i) + h * (k1 + 2 * k2 + 2 * k3 + k4) /6;
        fprintf('%10d%30.4f%30.4f%30.4f%30.4f%30.4f\n', i, y(1, i+1),y(2, i+1), t(i+1), dy(1), dy(2));
    end
end
