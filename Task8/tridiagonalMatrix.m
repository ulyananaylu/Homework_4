function [u] = tridiagonalMatrix(A, B, C, G)
    n = length(A) - 1;
    
    s = zeros(1, n + 1);
    t = zeros(1, n + 1);
    u = zeros(1, n + 1);
    
    s(1) = C(1) / B(1);
    t(1) = - G(1) / B(1);
    
    for i = 2:n + 1        
        s(i) = C(i) / (B(i) - A(i) * s(i-1));
        t(i) = (A(i) * t(i-1) - G(i)) / (B(i) - A(i) * s(i-1));
    end
    
    u(n + 1) = t(n + 1);
    
    for i = n:-1:1
        u(i) = s(i) * u(i + 1) + t(i);
    end
end
