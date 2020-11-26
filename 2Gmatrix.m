function [G] = Gmatrix(A)
    n = length(A);
    
    G = zeros(size(A));
    
    G(1,1) = sqrt(A(1,1));
    
    for i = 2:n
        G(1,i) = A(1,i) / G(1,1);
    end
    
    for i = 2:n
        sum = 0;
        for k = 1:i-1
            sum = sum + G(k,i) * G(k,i);
        end
        
        G(i,i) = sqrt(A(i,i) - sum);
        for j = i+1:n
            sum = 0;
            for k = 1:i-1
                sum = sum + G(k,i) * G(k,j);
            end
            
            G(i,j) = (A(i,j) - sum) / G(i,i);
        end
    end
    
    G = G';
    
end
