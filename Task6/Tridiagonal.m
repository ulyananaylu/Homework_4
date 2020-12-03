function [x, y] = Tridiagonal(a, b, h, alpha1, alpha2, alpha, beta1, beta2, beta)
    x = (a - h/2:h:b + h/2)';
    
    n = length(x) - 2;
    
    A = zeros(1, n + 2);
    B = zeros(1, n + 2);
    C = zeros(1, n + 2);
    G = zeros(1, n + 2);
    s = zeros(1, n + 2);
    t = zeros(1, n + 2);
    y = zeros(1, n + 2);
    
    B(1) = -alpha1/2 - alpha2/h;
    C(1) = alpha1/2 - alpha2/h;
    G(1) = alpha;
    
%     A(n + 2) = beta1/2 - beta2/h;
%     B(n + 2) = -beta1/2 - beta2/h;
%     G(n + 2) = beta;
    
    s(1) = C(1) / B(1);
    t(1) = - G(1) / B(1);
    
    for i = 2:n+1
        Pi = P(x(i));
        Qi = Q(x(i));
        Ri = R(x(i));
        Fi = F(x(i));
        
        A(i) = Pi / (h.^2) - Qi / (2 * h);
        B(i) = 2 * Pi / (h.^2) - Ri;
        C(i) = Qi / (2 * h) + Pi / (h.^2);
        G(i) = Fi;
        
        s(i) = C(i) / (B(i) - A(i) * s(i-1));
        t(i) = (A(i) * t(i-1) - G(i)) / (B(i) - A(i) * s(i-1));
    end
    
    y(n + 1) = t(n + 1);
    
    for i = n:-1:1
        y(i) = s(i) * y(i + 1) + t(i);
    end
end
