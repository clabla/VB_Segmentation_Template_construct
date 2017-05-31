%compute 1st and 2nd derivatives of lower bound---------------------------- 
%with respect to bias field parameters-------------------------------------

function [gb,Hb] = Bias_derivs(d,X,R,mu,Sigma,Bf,or,varargin)

if isempty(varargin)
    In.K = size(mu,2);
else
    In.K = varargin{1};
end

Bx = spm_dctmtx(d(1),or(1));
By = spm_dctmtx(d(2),or(2));
Bz = spm_dctmtx(d(3),or(3));
M  = size(Bx,2)*size(By,2)*size(Bz,2);
gb = zeros(M,1);
Hb = zeros(M,M);

for z = 1:d(3)
    
    i1  = d(1)*d(2)*(z-1) + 1;
    i2  = d(1)*d(2)*(z);
    
    x           = X(i1:i2);
    bf          = Bf(i1:i2);
    r           = R(i1:i2,:);
    r(isnan(r)) = 0;
    
    t   = zeros(numel(x),1);
    t1  = zeros(numel(x),1);
    
    i   = ~isnan(x) & x~=0;
        
    for k=1:In.K
        t(i)  = t(i) + bsxfun(@minus,x(i).*bf(i),mu(k)).*r(i,k)/Sigma(k);
        t1(i) = t1(i) + (x(i).*bf(i)).^2.*r(i,k)/Sigma(k);
    end
    
    t(i)  = t(i).*x(i).*bf(i) - 1;
    t1(i) = t1(i) + 1;
    
    Hb            = Hb + kron(Bz(z,:)'*Bz(z,:),spm_krutil(reshape(t1,d(1:2)),Bx,By,1)); 
    Hb(isnan(Hb)) = 0;
    gb            = gb + kron(Bz(z,:)',spm_krutil(reshape(t,d(1:2)),Bx,By,0)); 
    gb(isnan(gb)) = 0;
end
end
