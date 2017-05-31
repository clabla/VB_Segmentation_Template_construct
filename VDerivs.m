%compute 1st and 2nd derivatives of log likelihood-------------------------
%with respect to v                                              

function [gd,Hd] = VDerivs(dim,r,b,w,iM,ph,Jph,varargin)

if isempty(varargin)
    K   = size(r,2);
    prs = [1 1 1 0 0 0];
else
    K   = varargin{1};
    prs = varargin{2};
end
w   = reshape(w,[1,K]);
r   = single(reshape(r,[dim,K]));
gd  = zeros(prod(dim),3,'single');
db  = zeros(prod(dim),3,'single');
Hd  = zeros(prod(dim),6,'single');
E   = aff(ph,iM);
Jph = reshape(Jph,[prod(dim),3,3]);
bws = sum(bsxfun(@times,reshape(warp_template(b,[],E,prs),[prod(dim),K]),w),2);

for cl=1:K
    
    bt              = spm_diffeo('bsplins',single(log(b(:,:,:,cl))),E,prs);
    [~,dbx,dby,dbz] = spm_diffeo('bsplins',single(bt),Identity(dim,[1 1 1]),[2 2 2 0 0 0]);
    rt              = r(:,:,:,cl);
     
    j     = ~(~isfinite(bt) | ~isfinite(dbx)| ~isfinite(dby)| ~isfinite(dbz)| ~isfinite(rt) | (dbx==0 & dby==0 & dbz==0));
    
    bt    = w(cl)*exp(bt(j))./bws(j);
    

    dbx1 = [dbx(j),dby(j),dbz(j)]*iM(1:3,1).*Jph(j,1,1) + [dbx(j),dby(j),dbz(j)]*iM(1:3,2).*Jph(j,2,1) + [dbx(j),dby(j),dbz(j)]*iM(1:3,3).*Jph(j,3,1);
    dby1 = [dbx(j),dby(j),dbz(j)]*iM(1:3,1).*Jph(j,1,2) + [dbx(j),dby(j),dbz(j)]*iM(1:3,2).*Jph(j,2,2) + [dbx(j),dby(j),dbz(j)]*iM(1:3,3).*Jph(j,3,2);
    dbz1 = [dbx(j),dby(j),dbz(j)]*iM(1:3,1).*Jph(j,1,3) + [dbx(j),dby(j),dbz(j)]*iM(1:3,2).*Jph(j,2,3) + [dbx(j),dby(j),dbz(j)]*iM(1:3,3).*Jph(j,3,3);
    clear  dbx dby dbz
    
    gd(j,:) = gd(j,:) - bsxfun(@times,[dbx1,dby1,dbz1],(rt(j) - bt));
    
end

Hd(:,1) = Hd(:,1) + gd(:,1).*gd(:,1);
Hd(:,2) = Hd(:,2) + gd(:,2).*gd(:,2);
Hd(:,3) = Hd(:,3) + gd(:,3).*gd(:,3);
Hd(:,4) = Hd(:,4) + gd(:,1).*gd(:,2);
Hd(:,5) = Hd(:,5) + gd(:,1).*gd(:,3);
Hd(:,6) = Hd(:,6) + gd(:,2).*gd(:,3);

gd            = reshape(gd,[dim,3]);
gd(isnan(gd)) = 0;
Hd            = reshape(Hd,[dim,6]);
Hd(isnan(Hd)) = 0;
end
