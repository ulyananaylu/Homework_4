clear; clc; close('all');

%% Inputs:
%% Matrix
A = [
        -0.82005    -0.13542    0.26948;
        -0.13542    0.51486     0.002;
        0.26948     0.002     -0.83365
     ];

b = [
        1;
        2;
        3
     ];
 
 x = SimpleIteration(A, b, 0.00001);
 
 disp('Answer:');
 disp(x);
 
 x = Seidel(A, b, 0.00001);
 
 disp('Answer:');
 disp(x);

 
 disp('Answer by matlab:');
 disp(A \ b);
