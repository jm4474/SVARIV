function [hd,hdistanceMSW,hdistanceDMethod] = Hausdorff_Distance(n, horizons, MSWlbound, MSWubound,MSWlboundcum, MSWuboundcum, Dmethodlbound, Dmethodubound, Dmethodlboundcum, Dmethoduboundcum)
%-Reports the Hausdorff distance between MSW & Delta-Method confidence intervals
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
%       hd: Hausdorff distances computed for each variable and horizon (maximum of hdistanceMSW & hdistanceDMethod)  
%       hdistanceMSW: Distance starting from furthest point in MSW CI to closest point in DMethod CI 
%       hdistanceDMethod: Distance starting from furthest point in DMethod CI to closest point in MSW CI

%% Hausdorff Distance 

hdistanceMSW = Inf(n, horizons+1, 2); %taking x in MSW
hdistanceDMethod = Inf(n, horizons+1, 2); %taking x in D-method
hd = Inf(n, horizons+1, 2);

%Non-cumulative 

%For purposes of calculating distance from MSW to DMethod
a = MSWlbound;
b = MSWubound;
c = Dmethodlbound; 
d = Dmethodubound;

%For purposes of calculating distance from DMethod to MSW
e = Dmethodlbound;
f = Dmethodubound;
g = MSWlbound; 
h = MSWubound;

%Cumulative

%For purposes of calculating distance from MSW to DMethod
acum = MSWlboundcum;
bcum = MSWuboundcum;
ccum = Dmethodlboundcum; 
dcum = Dmethoduboundcum;

%For purposes of calculating distance from DMethod to MSW
ecum = Dmethodlboundcum;
fcum = Dmethoduboundcum;
gcum = MSWlboundcum; 
hcum = MSWuboundcum;

for var = 1:n
    
    for hori = 1:horizons+1
        %Starting in MSW CI, finding distance to closest point in DMethod CI
        
        %No intersections
        if b(var,hori) < c(var,hori) 

            hdistanceMSW(var,hori,1) = abs(c(var,hori) - a(var,hori)); 

        end

        if a(var,hori) > d(var,hori) 

            hdistanceMSW(var,hori,1) = abs(b(var,hori) - d(var,hori)); 

        end
        
        %Intersections
        if b(var,hori) >= c(var,hori) && b(var,hori) <= d(var,hori) && a(var,hori) < c(var,hori)   

            hdistanceMSW(var,hori,1) = abs(a(var,hori) - c(var,hori)); 

        end 

        if a(var,hori) >= c(var,hori) && a(var,hori) <= d(var,hori) && b(var,hori) > d(var,hori) 

            hdistanceMSW(var,hori,1) = abs(b(var,hori) - d(var,hori));

        end 

        if a(var,hori) >= c(var,hori) && a(var,hori) <= d(var,hori) && b(var,hori) >= c(var,hori) && b(var,hori) <= d(var,hori) 

            hdistanceMSW(var,hori,1) = 0;

        end

        if a(var,hori) <= c(var,hori) && b(var,hori) >= d(var,hori) 

            hdistanceMSW(var,hori,1) = max(abs(c(var,hori)-a(var,hori)),abs(d(var,hori)-b(var,hori))); 

        end
        
        %Starting in DMethod CI, finding distance to closest point in MSW CI

        if f(var,hori) < g(var,hori) %Case where [a,b] (MSW CI) < [c,d](D-Method CI) (no intersection)

            hdistanceDMethod(var,hori,1) = abs(g - e); 

        end

        if e(var,hori) > h(var,hori) %Case where [a,b] > [c,d] (no intersection)

            hdistanceDMethod(var,hori,1) = abs(f(var,hori) - h(var,hori)); 

        end

        if f(var,hori) >= g(var,hori) && f(var,hori) <= h(var,hori) && e(var,hori) < g(var,hori)  %Case where c <= b <= d 

            hdistanceDMethod(var,hori,1) = abs(e(var,hori) - g(var,hori)); 

        end 

        if e(var,hori) >= g(var,hori) && e(var,hori) <= h(var,hori) && f(var,hori) > h(var,hori) %Case where c <= a <= d 

            hdistanceDMethod(var,hori,1) = abs(f(var,hori) - h(var,hori));

        end 

        if e(var,hori) >= g(var,hori) && e(var,hori) <= h(var,hori) && f(var,hori) >= g(var,hori) && f(var,hori) <= h(var,hori) 

            hdistanceDMethod(var,hori,1) = 0;

        end

        if e(var,hori) <= g(var,hori) && f(var,hori) >= h(var,hori) 

            hdistanceDMethod(var,hori,1) = max(abs(g(var,hori)-e(var,hori)),abs(h(var,hori)-f(var,hori))); 

        end

        %Cumulative 
        
        %Starting in MSW CI, finding distance to closest point in DMethod CI
        
        %No intersections

        if bcum(var,hori) < ccum(var,hori)

            hdistanceMSW(var,hori,2) = abs(ccum(var,hori) - acum(var,hori)); 

        end

        if acum(var,hori) > dcum(var,hori) 

            hdistanceMSW(var,hori,2) = abs(bcum(var,hori) - dcum(var,hori)); 

        end
        
        %Intersections

        if bcum(var,hori) >= ccum(var,hori) && bcum(var,hori) <= dcum(var,hori) && acum(var,hori) < ccum(var,hori)   

            hdistanceMSW(var,hori,2) = abs(acum(var,hori) - ccum(var,hori)); 

        end 

        if acum(var,hori) >= ccum(var,hori) && acum(var,hori) <= dcum(var,hori) && bcum(var,hori) > dcum(var,hori)  

            hdistanceMSW(var,hori,2) = abs(bcum(var,hori) - dcum(var,hori));

        end 

        if acum(var,hori) >= ccum(var,hori) && acum(var,hori) <= dcum(var,hori) && bcum(var,hori) >= ccum(var,hori) && bcum(var,hori) <= dcum(var,hori) 

            hdistanceMSW(var,hori,2) = 0;

        end

        if acum(var,hori) <= ccum(var,hori) && bcum(var,hori) >= dcum(var,hori) 

            hdistanceMSW(var,hori,2) = max(abs(ccum(var,hori)-acum(var,hori)),abs(dcum(var,hori)-bcum(var,hori))); 

        end 

        %Starting in DMethod CI, finding distance to closest point in MSW CI 
        
        %No intersections
        
        if fcum(var,hori) < gcum(var,hori) 

            hdistanceDMethod(var,hori,2) = abs(gcum(var,hori) - ecum(var,hori)); 

        end

        if ecum(var,hori) > hcum(var,hori) 

            hdistanceDMethod(var,hori,2) = abs(fcum(var,hori) - hcum(var,hori)); 

        end
        
        %Intersections

        if fcum(var,hori) >= gcum(var,hori) && fcum(var,hori) <= hcum(var,hori) && ecum(var,hori) < gcum(var,hori)   

            hdistanceDMethod(var,hori,2) = abs(ecum(var,hori) - gcum(var,hori)); 

        end 

        if ecum(var,hori) >= gcum(var,hori) && ecum(var,hori) <= hcum(var,hori) && fcum(var,hori) > hcum(var,hori) 

            hdistanceDMethod(var,hori,2) = abs(fcum(var,hori) - hcum(var,hori));

        end 

        if ecum(var,hori) >= gcum(var,hori) && ecum(var,hori) <= hcum(var,hori) && fcum(var,hori) >= gcum(var,hori) && fcum(var,hori) <= hcum(var,hori) 

            hdistanceDMethod(var,hori,2) = 0;

        end

        if ecum(var,hori) <= gcum(var,hori) && fcum(var,hori) >= hcum(var,hori) 

            hdistanceDMethod(var,hori,2) = max(abs(gcum(var,hori)-ecum(var,hori)),abs(hcum(var,hori)-fcum(var,hori))); 

        end
        
    end
    
end

hd(:,:,:) = max(hdistanceMSW,hdistanceDMethod);
    
end 
