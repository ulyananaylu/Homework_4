function [t, y] = BackwardEuler(t0, y0, h, tn, A);
    W = inv(eye(size(A))-A*h);
    t = (t0:h:tn)';
    s = size(A);
    y = zeros(s(1), length(t));
    
    y(:,1) = y0;
    for i = 1:1:length(t)-1
        y(:,i+1) = W * y(:,i);
        
    end
end
