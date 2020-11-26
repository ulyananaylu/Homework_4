clear; clc; close('all');

%% Inputs:
%% Matrix
A = [
        -0.82005    -0.13542    0.26948;
        -0.13542    0.51486     0.2;
        0.26948     0.2     -0.83365
     ];
 
 [X, lambda] = Jacobi(A, 0.0000001);
 format long;
 disp('Eigenvectors by Jacobi:');
 disp(X);
 disp('Eigenvalues by Jacobi:');
 disp(lambda);
 
 [R,D] = eig(A); 
 disp('Eigenvectors by matlab:');
 disp(R); 
 disp('Eigenvalues by matlab:');
 disp(diag(D)');
