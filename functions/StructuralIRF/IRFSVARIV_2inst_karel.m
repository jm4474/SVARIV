function IRFs = IRFSVARIV_2inst_karel(AL,Sigma,Gamma,hori,x,nvar)

% Note k is the number of shocks with instruments available

n  = length(Sigma); 

p  = size(AL,2)/n;

k  = size(Gamma,2);
 
Gammap1     = Gamma(1:k,1:k)';     % Sigmamu1' in AER notation

Gammap2     = Gamma(k+1:n,1:k)';   % Sigmamu2' 

b21ib11     = (Gammap1\Gammap2)';  % Sigmamu2*(Sigmamu1^{-1}) 

Sig11       = Sigma(1:k,1:k);      

Sig21       = Sigma(k+1:n,1:k); 

Sig22       = Sigma(k+1:n,k+1:n);

ZZp         = b21ib11*Sig11*b21ib11'-(Sig21*b21ib11'+b21ib11*Sig21')+Sig22;

b12b12p     = (Sig21- b21ib11*Sig11)'*(ZZp\(Sig21- b21ib11*Sig11));
b11b11p     = Sig11-b12b12p;
b22b22p     = Sig22+b21ib11*(b12b12p-Sig11)*b21ib11';
b12ib22     = ((Sig21- b21ib11*Sig11)'+b12b12p*b21ib11')/(b22b22p');
b11iSig     = eye(k)/(eye(k)-b12ib22*b21ib11);
b21iSig     = b21ib11*b11iSig;

SigmaTSigmaTp = b11iSig\b11b11p/b11iSig';

if min(eig(Sigma))>0
    
% This is neccessary because otherwise we get a warning message

SigmaT = chol(SigmaTSigmaTp,'lower');

DT     = [b11iSig;b21iSig]*SigmaT;

DT     = x*DT(:,1)./DT(nvar,1);

Cauxsim= [eye(n),MARep(AL,p,hori)]; 

Csim   = reshape(Cauxsim,[n,n,hori+1]);

IRFs = reshape(sum(bsxfun(@times,Csim,DT'),2),[n,hori+1]);

% The implementation based on Cauxsim and Csim can be replaced by:
% irs = [];
%  irs(p+1,:) = DT;
%  for tt=2:hori+1
%  lvars = (irs(p+tt-1:-1:tt,:))';
%  irs(p+tt,:) = lvars(:)'*AL(:,1:p*n)';     
%  end
% 
% IRFs = irs(p+1:end,:)'; 

else 
    
IRFs = NaN(n,hori+1); 

end

end
