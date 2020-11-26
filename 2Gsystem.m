function [x] = Gsystem(A, b)
    G = Gmatrix(A);
    
    fprintf('GSystem\n\n');
    
    n = size(b);
    
    y = zeros(n);
    
    y(1) = b(1) / G(1,1);
    for i=2:n
        sum = 0;
        for k = 1:i-1
            sum = sum + G(i,k)*y(k);
        end
        
        y(i) = (b(i) - sum)/G(i,i);
    end
    
    x = zeros(n);
    
    n = length(b);
    
    x(n) = y(n) / G(n,n);
    
    G = G';
    
    for i = n-1:-1:1
        sum = 0;
        for k = i+1:n
            sum = sum + G(i,k)*x(k);
        end
        x(i) = (y(i) - sum) / G(i,i);
    end
end
