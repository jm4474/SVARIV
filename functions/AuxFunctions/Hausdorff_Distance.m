function [hd,hdistance,hdistance2] = Hausdorff_Distance(n, horizons, MSWlbound, MSWubound,MSWlboundcum, MSWuboundcum, Dmethodlbound, Dmethodubound, Dmethodlboundcum, Dmethoduboundcum)
%-Reports the Hausdorff distance between MSW & Delta-Method confidence
%intervals
%    -Syntax:
%       [hd,hdistance,hdistance2] = Hausdorff_Distance(n, horizons, MSWlbound, MSWubound,MSWlboundcum, MSWuboundcum, Dmethodlbound, Dmethodubound, Dmethodlboundcum, Dmethoduboundcum)
%    -Inputs:
%       n:                   Number of variables in the VAR model 
%       horizons:            Number of horizons for the Impulse Response Functions (IRFs)               (1 x 1)  
%                            (does not include the impact horizon 0)    
%       MSWlbound:           lower bounds of MSW confidence interval                                    (n x horizons+1)
%       MSWubound:           upper bounds of MSW confidence interval                                    (n x horizons+1)
%       MSWlboundcum:        lower bounds of MSW cumulative confidence interval                         (n x horizons+1)
%       MSWuboundcum:        upper bounds of MSW cumulative confidence interval                         (n x horizons+1)
%       Dmethodlbound:       lower bounds of delta method cumulative confidence interval                (n x horizons+1)
%       Dmethodubound:       upper bounds of delta method cumulative confidence interval                (n x horizons+1)
%       Dmethodlboundcum:    lower bounds of delta method cumulative confidence interval                (n x horizons+1)
%       Dmethoduboundcum:    upper bounds of delta method cumulative confidence interval                (n x horizons+1)
%
%   -Outputs:
%       hd: Hausdorff distances computed for each variable (this was to
%           check the way Pepe suggested doing it)
%       hdistance: Grid of Hausdorff distances for each variable and each
%           horizon (likely not the correct way, keeping just in case) 
%       hdistance2: Hausdorff distances computed for each variable (I think
%           the way Pepe suggested doing it)


%% Initial way of trying it out, creating Hausdorff distance for each
% %variable, each horizon (THIS IS PROBABLY WRONG!!!)

hdistance = zeros(n, horizons+1, 2);

hdistance(:,:,1) = max(abs(MSWlbound-Dmethodlbound), abs(MSWubound-Dmethodubound)); 

hdistance(:,:,2) = max(abs(MSWlboundcum-Dmethodlboundcum), abs(MSWuboundcum-Dmethoduboundcum));

%% Second way of trying it out, creating one Hausdorff distance for each
%variable

hdistance2 = zeros(n,2); 

for var = 1:n 
    
    %Calculating max and min bounds for non-cumulative MSW & Dmethod for each variable and
    %then taking hdistance
    
    MSWlbound2 = min(MSWlbound(var,:)); 
    
    MSWubound2 = max(MSWubound(var,:)); 
    
    Dmethodlbound2 = min(Dmethodlbound(var,:));
    
    Dmethodubound2 = max(Dmethodubound(var,:));
    
    hdistance2(var,1) = max(abs(MSWlbound2 - Dmethodlbound2), abs(MSWubound2 - Dmethodubound2));
    
    
    %Calculating max and min bounds for cumulative MSW & Dmethod for each variable and
    %then taking hdistance
    MSWlboundcum2 = min(MSWlboundcum(var,:)); 
    
    MSWuboundcum2 = max(MSWuboundcum(var,:)); 
    
    Dmethodlboundcum2 = min(Dmethodlboundcum(var,:));
    
    Dmethoduboundcum2 = max(Dmethoduboundcum(var,:));
    
    hdistance2(var,2) = max(abs(MSWlboundcum2 - Dmethodlboundcum2), abs(MSWuboundcum2 - Dmethoduboundcum2)); 
    
end


%% 3rd way of trying it out 
%Non-cumulative 

MSWbounds = horzcat(MSWlbound, MSWubound); 

Dmethodbounds = horzcat(Dmethodlbound, Dmethodubound); 

%Distances from MSWbounds to Dmethodbounds
d = zeros(size(MSWbounds,2),size(Dmethodbounds,2),n); 

shortest = Inf(size(MSWbounds,2),n);

h1 = zeros(n,1);

for var = 1:n 
    
    for i = 1:size(MSWbounds,2) 
        
        for j = 1:size(Dmethodbounds,2)
            
            d(i,j,var) = abs(MSWbounds(var,i) - Dmethodbounds(var,j));
            
        end   
        
        if min(d(i,:,var)) < shortest(i,var)

            shortest(i,var) = min(d(i,:,var));

        else           
        end
 
    end

    h1(var,1) = max(shortest(:,var)); 
    
end

%Distances from Dmethodbounds to MSWbounds
d = zeros(size(MSWbounds,2),size(Dmethodbounds,2),n); 

shortest = Inf(size(MSWbounds,2),n);

h2 = zeros(n,1);

for var = 1:n 
    
    for i = 1:size(Dmethodbounds,2) 
        
        for j = 1:size(MSWbounds,2)
            
            d(i,j,var) = abs(Dmethodbounds(var,i) - MSWbounds(var,j));
            
        end   
        
        if min(d(i,:,var)) < shortest(i,var)

            shortest(i,var) = min(d(i,:,var));

        else           
        end
 
    end

    h2(var,1) = max(shortest(:,var)); 
    
end

%Hausdorff distances BETWEEN bounds (non-cumulative)
hd(:,1) = bsxfun(@max,h1,h2);    

%Cumulative 

MSWboundscum = horzcat(MSWlboundcum, MSWuboundcum); 

Dmethodboundscum = horzcat(Dmethodlboundcum, Dmethoduboundcum); 

d = zeros(size(MSWboundscum,2),size(Dmethodboundscum,2),n); 

shortest = Inf(size(MSWboundscum,2),n);

h3 = zeros(n,1);

for var = 1:n 
    
    for i = 1:size(MSWboundscum,2) 
        
        for j = 1:size(Dmethodboundscum,2)
            
            d(i,j,var) = abs(MSWboundscum(var,i) - Dmethodboundscum(var,j));
            
        end   
        
        if min(d(i,:,var)) < shortest(i,var)

            shortest(i,var) = min(d(i,:,var));

        else           
        end
 
    end

    h3(var,1) = max(shortest(:,var)); 
    
end

d = zeros(size(MSWboundscum,2),size(Dmethodboundscum,2),n); 

shortest = Inf(size(MSWboundscum,2),n);

h4 = zeros(n,1);

for var = 1:n 
    
    for i = 1:size(Dmethodboundscum,2) 
        
        for j = 1:size(MSWboundscum,2)
            
            d(i,j,var) = abs(Dmethodboundscum(var,i) - MSWboundscum(var,j));
            
        end   
        
        if min(d(i,:,var)) < shortest(i,var)

            shortest(i,var) = min(d(i,:,var));

        else           
        end
 
    end

    h4(var,1) = max(shortest(:,var)); 
    
end

%Hausdorff distances BETWEEN bounds (cumulative)
hd(:,2) = bsxfun(@max,h3,h4);    

    
end 
