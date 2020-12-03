function [L, U] = LU(A)
    L = zeros(size(A));
    U = zeros(size(A));
    
    n = length(A);
    
    for i = 1:n
        for j = i:n
            sumL = 0;
            sumU = 0;
            for k = 1:i-1
                sumL = sumL + L(j,k)*U(k,i);
                sumU = sumU + L(i,k)*U(k,j);
            end
            L(j,i) = A(j,i) - sumL;
            U(i,j) = (A(i,j) - sumU)/L(i,i);
        end
    end
end
