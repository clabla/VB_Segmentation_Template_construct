%compute responsabilities (VB)---------------------------------------------

function q = resp(x,bf,mp,po,varargin)

C  = size(x,2);
K  = size(mp,2);
q  = NaN(size(x,1),K);
li = bi2de(double(~isnan(x) & x~=0)) + 1;
mp(isnan(mp)) = 1/K;
mp(mp<=0)      = exp(-1e+1);

if ~isempty(varargin)
    K12 = varargin{1};
else
    K12 = 1:K;
end

for p=2:2^C
    
    ci = logical(de2bi(p-1,C));
    ri = repmat(p,[size(x,1),1])==li;
    D1  = sum(ci);
    
    for k=K12
        ldW     = 2*sum(log(diag(chol(po.W(ci,ci,k)))));
        Sdig    = sum(psi((1 + po.nu(1,k)-(1:D1))/2));
        d       = bsxfun(@minus,x(ri,ci).*bf(ri,ci),po.m(ci,k)')*chol(po.W(ci,ci,k));
        EQ      = D1/po.beta(1,k) + po.nu(1,k)*sum(d.*d,2);
        q(ri,k) = log(prod(bf(ri,ci),2)) + log(mp(ri,k)) - D1/2*log(2*pi) + 0.5*(ldW + D1*log(2) + Sdig) - 0.5*EQ;
    end
end

q = bsxfun(@minus,q,nanmax(q,[],2));
q = bsxfun(@rdivide,exp(q),nansum(exp(q),2));
end

