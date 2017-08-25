function [M] = vechtovec(n)
%Matrix M such that vec(Sigma) = M*vech(Sigma) where Sigma is n x n

%% 1) Initialize values

M = zeros(n,(n)*(n+1)/2);

for m = 1:n
   
    auxcol1 = ((m-1)*(n+1)) - ((m-1)*(m)/2); 
    % sum from j=1 to m-1 of all blocks of size n+1-j
    
    auxcol2 = ((m)*(n+1)) - ((m)*(m+1)/2);
    % sum from j=1 to m of all blocks of size n+1-j
    
    aux     = eye(n+1-m);
    
    M(n*(m-1)+(m-1)+1:n*m, auxcol1+1:auxcol2 ) = aux;
    
    for j = 1:n-m
        
    M(((m+j-1)*n)+m, auxcol1+1:auxcol2 )  = aux(j+1,:);   
    
    end
  
   
end

end

