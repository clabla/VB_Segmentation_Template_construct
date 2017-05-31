%compute 1st and 2nd derivatives of lower bound----------------------------
%with respect to affine deformation parameters (using gradients of b)

function [gd,Hd] = AffDerivs(r,b,iM,diM,ph,C,a,w,varargin)

if isempty(varargin)
    K   = size(r,2);
    prs = [1 1 1 0 0 0];
else
    K   = varargin{1};
    prs = varargin{2};
end

gd  = zeros(1,12);
Hd  = zeros(12,12);
d   = size(ph(:,:,:,1));
r   = single(reshape(r,[d,K]));
E   = aff(ph,iM);
bws = sum(bsxfun(@times,reshape(warp_template(b,[],E,prs),[prod(d),K]),w),2);
x   = E(:,:,:,1);
y   = E(:,:,:,2);
z   = E(:,:,:,3);
clear E

for cl=1:K
    
     bt              = spm_diffeo('bsplins',single(log(b(:,:,:,cl))),cat(4,x,y,z),prs);
     [~,dbx,dby,dbz] = spm_diffeo('bsplins',single(bt),Identity(d,[1 1 1]),[2 2 2 0 0 0]);
     rt              = r(:,:,:,cl);
     
     %exclude voxels where tpms are not defined or where gradients are zero
     j  = ~(~isfinite(bt) | ~isfinite(dbx)| ~isfinite(dby)| ~isfinite(dbz)| ~isfinite(rt) | (dbx==0 & dby==0 & dbz==0));
       
     bt = w(cl)*exp(bt(j))./bws(j);    
     
     A  =  [x(j).*dbx(j) y(j).*dbx(j) z(j).*dbx(j) dbx(j) ...
            x(j).*dby(j) y(j).*dby(j) z(j).*dby(j) dby(j) ...
            x(j).*dbz(j) y(j).*dbz(j) z(j).*dbz(j) dbz(j)];             
     A  = (A*diM);
     
     gd = gd - sum(bsxfun(@times,A,(rt(j) - bt)),1);    
     Hd = Hd + bsxfun(@times,A,(rt(j) - bt))'*bsxfun(@times,A,(rt(j) - bt));
end

gd = gd'+ C\a;
Hd = Hd + inv(C);
gd = double(gd);
Hd = double(Hd);
end


