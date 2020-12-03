function [t, y] = Adams3(t0, y0, h, tn, A)
    W1 = inv((eye(size(A))-A*h*5/12)) * (eye(size(A))+2*h*A/3);
    W2 = inv((eye(size(A))-A*h*5/12)) * (h*A/12);
    fprintf('Third order Adams method\n\n');
    fprintf('Eigenvalues of A matrix: %1.4f and %1.4f\n\n', eig(A).');
    fprintf('%10s%30s%30s%30s%30s%30s\n', 'i', 'y1i', 'y2i', 'ti', 'dy1/dt(yi,ti)', 'dy2/dt(yi,ti)');
    dy = A * y0;
    fprintf('%10d%30.4f%30.4f%30.4f%30.4f%30.4f\n', 0, y0(1), y0(2), t0, dy(1), dy(2));
    t = (t0:h:tn)';
    y = zeros(2, length(t));
    
    y(:,1) = y0;
    
    [tt, backwardEuler] = BackwardEuler(t0, y0, h, t0+h, A, false);
    
    y(:,2) = backwardEuler(:,2);
    
    fprintf('%10d%30.4f%30.4f%30.4f%30.4f%30.4f\n', 1, y(1, 2),y(2, 2), t(2), dy(1), dy(2));
    
    for i = 2:length(t)-1
        y(:,i+1) = W1 * y(:,i) - W2 * y(:, i-1);
        fprintf('%10d%30.4f%30.4f%30.4f%30.4f%30.4f\n', i, y(1, i+1),y(2, i+1), t(i+1), dy(1), dy(2));
    end
end
