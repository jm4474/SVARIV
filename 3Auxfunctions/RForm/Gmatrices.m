function [G,Gcum] = Gmatrices(AL,C,p,hori,n)
% -------------------------------------------------------------------------
% Computes the derivatives of vec(C) wrt vec(A) based on
% Lütkepohl H. New introduction to multiple time series analysis. ? Springer, 2007.
% 
% Inputs:
% - AL: VAR model coefficients
% - C: MA representation coefficients
% - p: lag order
% - hori: forecast horizon
% - n: number of variables
% Outputs:
% - G: derivatives of C wrt A 
% 
% This version: March 31, 2015
% Please, cite Gafarov, B. and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
% -------------------------------------------------------------------------


%% A and J matrices in Lutkepohl's formula for the derivative of C with respect to A
J = [eye(n), zeros(n,(p-1)*n)];
Alut = [AL; eye(n*(p-1)),zeros(n*(p-1),n)];


%% AJ is a 3D array that contains A^(k-1) J' in the kth 2D page of the the 3D array
AJ = zeros(n*p, n, hori);
for k=1:hori
    AJ(:,:,k) = ((Alut)^(k-1)) * J';
end


%% matrix [ JA'^0; JA'^1; ... J'A^{k-1} ];
JAp = reshape(AJ, [n*p,n*hori])';


%% G matrices
AJaux = zeros(size(JAp,1)*n, size(JAp,2)*n, hori);
Caux = reshape([eye(n), C(:,1:(hori-1)*n)], [n,n,hori]);
for i=1:hori
    AJaux(((n^2)*(i-1))+1:end,:,i) = kron(JAp(1:n*(hori+1-i),:), Caux(:,:,i)); 
end
Gaux = permute(reshape(sum(AJaux,3)', [(n^2)*p, n^2, hori]), [2,1,3]);
G = zeros(size(Gaux,1), size(Gaux,2), size(Gaux,3)+1);
G(:,:,2:end) = Gaux;


%% Cumulative version of G matrices
Gcum = cumsum(G,3);


end
