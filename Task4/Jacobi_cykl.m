function lambda = Jacobi_cykl(A, eps)
    AbsA = abs(triu(A,1));
    

    n = length(A);

    
    lambda = zeros(1, n);
    
    num = 0;
    m = max(max(AbsA));
    while m > eps
        for Ik = 1:n
            for Jk = Ik+1:n
            Ak = zeros(size(A));
            
            d = sqrt((A(Ik,Ik) - A(Jk,Jk)).^2 +4*A(Ik,Jk).^2);
            c = sqrt((1+abs(A(Ik,Ik) - A(Jk,Jk))/d)/2);
            s = sign(A(Ik,Jk) * (A(Ik,Ik) - A(Jk,Jk))) * sqrt((1 - abs(A(Ik,Ik) - A(Jk,Jk))/d)/2);
        
            for i = 1:n
                for j = 1:n
                    if (i ~= Ik) && (i ~= Jk) && (j ~= Ik) && (j ~= Jk)
                        Ak(i,j) = A(i,j);
                        
                    end
                
                    if (i ~= Ik) && (i ~= Jk) 
                        Ak(i,Ik) = c * A(i,Ik) + s * A(i,Jk);
                        Ak(Ik,i) = c * A(i,Ik) + s * A(i,Jk);
                    
                        Ak(i,Jk) = c * A(i,Jk) - s * A(i,Ik);
                        Ak(Jk,i) = c * A(i,Jk) - s * A(i,Ik);
                    
                        
                    end
                
                end
            end
        
            Ak(Ik, Ik) = A(Ik,Ik) * (c.^2) + 2*c*s*A(Ik,Jk) + (s.^2) * A(Jk,Jk);
            Ak(Jk, Jk) = (s.^2) * A(Ik,Ik) - 2*c*s*A(Ik,Jk) + (c.^2) * A(Jk,Jk);
                    
            Ak(Ik,Jk) = 0;
            Ak(Jk,Ik) = 0;        
        

             A = Ak;
        
            AbsA = abs(triu(A,1));
            m = max(max(AbsA));

        
            num = num + 1;
            end
        end
    end
    disp('Количество итераций циклическим методом: ')
    disp(num);
    for i = 1:n
        sum = 0;
        for j = i:n
            if i ~= j
                sum = sum + A(i,j)/(A(i,i) - A(j,j));
            end
        end
        
        lambda(i) = A(i,i) + sum;
    end
end
