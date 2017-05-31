function [po,pr,mom,Lk,Elpz,Eqz,r,w]=VB_EM(Xm,Bf,bs,po,pr,w,K,N,in,varargin)

if ~isempty(varargin)
    rl    = varargin{1};   
    K1    = size(rl,2);  
    j     = rl<0.1;
    rl(j) = 0;
    rl    = bsxfun(@rdivide,rl,nansum(rl,2)+eps);
end
bs = reshape(bs,[N,K]);
r  = zeros(N,K);
s0 = zeros(1,K);
if isempty(po)
    r = bs;
else
    if isempty(varargin)
        r = Resp(Xm,Bf,bs,po);
    else
        if isempty(varargin{1})
            r = Resp(Xm,Bf,bs,po);
        else
            if isempty(varargin{2})
                t = 1:K1;
            else
                t  = varargin{2};
            end
            j = sum(rl,2)<=0.5;
            for kt = t
                k1 = in == kt;
                k2 = t == kt;
                i      = rl(:,k2)>0.1;
                if any(i)
                    r(i,:) = nansum(cat(3,r(i,:),Resp(Xm(i,:),Bf(i,:),bs(i,:),po,k1)),3);
                end
            end
            r(j,:)  = nansum(cat(3,r(j,:),Resp(Xm(j,:),Bf(j,:),bs(j,:),po)),3);
        end
    end
end
r(r==0) = NaN;
r       = bsxfun(@minus,r,nanmax(r,[],2));
r       = bsxfun(@rdivide,exp(r),nansum(exp(r),2));
Elpz    = 0;
Eqz     = 0;
mom     = spm_SuffStats(Xm.*Bf,r);
for k=1:K,
    Elpz = Elpz + nansum(log(bs(:,k)+eps).*r(:,k));
    Eqz  = Eqz  + nansum(log(r(:,k)+eps).*r(:,k));
end
[po,Lk,pr] = spm_VBGaussiansFromSuffStats(mom,pr);

for m = 1:numel(mom)
    s0 = s0 + mom(m).s0;
end
for k=1:K
    w(k) = (s0(k))/(nansum(bs(:,k)) + eps);
end
end

%compute responsabilities (VB)---------------------------------------------

function q = Resp(x,bf,mp,po,varargin)

C  = size(x,2);
K  = size(mp,2);
q  = NaN(size(x,1),K);
li = bi2de(double(~isnan(x) & x~=0)) + 1;
mp(mp<=0) = eps;

if ~isempty(varargin)
    K12 = find(varargin{1});
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
end

