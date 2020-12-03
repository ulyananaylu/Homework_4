function tmp = myjacobi(k, n)
    pj=sym(zeros(n+1,1));
    syms x;
    pj(1)=1;
    
    if n == 1
        pj(2) = (k + 1) * x;
    elseif n > 1
        pj(2) = (k + 1) * x;
        
        for i = 2:n
            pj(i+1)=((i+k)/(i+2*k))*(2+(2*k-1)/i)*x*pj(i)-...
                ((i+k)/(i+2*k))*(1+(k-1)/i)*pj(i-1);
        end
    end
   
    tmp = pj;
end
