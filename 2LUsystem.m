function [x] = LUsystem(A,b)
    [L, U] = LU(A);
    fprintf('LUSystem\n\n');
    n = size(b);
    
    y = zeros(n);
    
    y(1) = b(1) / L(1,1);
    for i=2:n
        sum = 0;
        for k = 1:i-1
            sum = sum + L(i,k)*y(k);
        end
        
        y(i) = (b(i) - sum)/L(i,i);
    end
    
    x = zeros(n);
    
    n = length(b);
    
    x(n) = y(n);
    
    for i = n-1:-1:1
        sum = 0;
        for k = i+1:n
            sum = sum + U(i,k)*x(k);
        end
        x(i) = y(i) - sum;
    end
    
end
