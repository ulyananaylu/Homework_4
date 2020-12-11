function [t, y] = Adams3(t0, y0, h, tn, A);
    W1 = inv((eye(size(A))-A*h*5/12)) * (eye(size(A))+2*h*A/3);
    W2 = inv((eye(size(A))-A*h*5/12)) * (h*A/12);
    t = (t0:h:tn)';
    s = size(A);
    y = zeros(s(1), length(t));
    
    y(:,1) = y0;
    
    [~, backwardEuler] = BackwardEuler(t0, y0, h, t0+h, A);
    
    y(:,2) = backwardEuler(:,2);
    
    for i = 2:length(t)-1
        y(:,i+1) = W1 * y(:,i) - W2 * y(:, i-1);
    end
end
