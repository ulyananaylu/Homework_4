function [t, y] = Euler(t0, y0, h, tn, A)
    W = eye(size(A))+A*h;
    fprintf('Euler method\n\n');
    fprintf('Eigenvalues of W matrix: %1.4f and %1.4f\n\n', eig(W).');
    fprintf('%10s%30s%30s%30s%30s%30s\n', 'i', 'y1i', 'y2i', 'ti', 'dy1/dt(yi,ti)', 'dy2/dt(yi,ti)');
    dy = A * y0;
    fprintf('%10d%30.4f%30.4f%30.4f%30.4f%30.4f\n', 0, y0(1), y0(2), t0, dy(1), dy(2));
    t = (t0:h:tn)';
    y = zeros(2, length(t));
    
    y(:,1) = y0;
    for i = 1:1:length(t)-1
        y(:,i+1) = W * y(:,i);
        fprintf('%10d%30.4f%30.4f%30.4f%30.4f%30.4f\n', i, y(1, i+1),y(2, i+1), t(i+1), dy(1), dy(2));
    end
end
