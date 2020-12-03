A=[
-0.951210  -0.097787   0.358427
-0.097787   0.615454    0.022286
0.358427    0.022286    -0.957290];

Yn=[1;1;1;];
l1=0;
l2=1e100;
j=0;
while abs(l1-l2)>0.001
    Y=Yn;
    l2=l1;
    Yn=A*Y;
    Yn'
    l1=Yn(1)/Y(1)
    j=j+1;
    l(j)=l1;
end
Y
j
eigs(A,1)%самое большое

j=2;
la(1)=1e100;la(2)=0;
while abs(la(j)-la(j-1))>0.001
j=j+1;
la(j)=(l(j)*l(j+2)-(l(j+1))^2)/(l(j)-2*l(j+1)+l(j+2));
end
j
la(j)

j=2; 
Yn=[1;1;1]; 
la(1)=1e100;la(2)=0; 
while abs(la(j)-la(j-1))>0.001 
Y=A*Yn; 
j=j+1; 
la(j)=(Y'*Yn)/(Yn'*Yn); 
Yn=Y; 
end 
j 
la(j)
