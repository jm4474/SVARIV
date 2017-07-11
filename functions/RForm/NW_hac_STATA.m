function Sigma= NW_hac_STATA(vars,lags)

Sigma0=(1/(size(vars,1)))*(vars'*vars);

Sigma_cov=@(k) (1/(size(vars,1)))...
               *(vars(1:end-k,:))'...
               *(vars(1+k:end,:));

Sigma=Sigma0;
for n=1:lags
    Sigma=Sigma+(1-n/(lags+1))*(Sigma_cov(n)+Sigma_cov(n)');
end
end