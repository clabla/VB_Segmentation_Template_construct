function [fntu, msk] = Init_u(DG,MG,C,K1,S,fn,fnr,udir,vdir,mm,affpar,prm)

fntu   = cell(S,1);
Ker    = spm_shoot_greens('kernel',DG,prm);
msk    = zeros([DG,K1]);

parfor s=1:S
    u    = [];
    tmp  = zeros([DG,K1]);
    if ~isempty(udir)
        for c=1:C
            [~,na] = fileparts(fn{s,c});
            if exist(fullfile(udir,strcat('u',na,'.nii')),'file')
                U = nifti(fullfile(udir,strcat('u',na,'.nii')));
                u = warp(single(U.dat(:,:,:,:)),U.mat\MG,Identity(DG,[1 1 1]),[2 2 2 1 1 1]);
%                 u = bsxfun(@times,u,reshape(sqrt(sum(U.mat(1:3,1:3).^2))./mm,[1 1 1 3]));
                u(isnan(u)) = 0;
                break
            else
                u = zeros([DG,3],'single');
                break
            end
        end
    else
        u = zeros([DG,3],'single');
    end
    
    if ~isempty(fnr{s})
        
        R             = nifti(fnr{s}.name);
        Kr            = fnr{s}.K;
        r             = single(R.dat(:,:,:,:));
        [~,~,~,th]    = spm_shoot3d(u,prm,6,Ker);
        th            = spm_diffeo('comp',th,aff(Identity(DG,[1 1 1]),Mtransf(affpar{s},MG,MG)));
        if numel(Kr)==size(r,4)
            tmp(:,:,:,Kr) = Warp_template(r,eye(4),th,[1 1 1 0 0 0]);
        else
            tmp(:,:,:,Kr) = repmat(Warp_template(r,eye(4),th,[1 1 1 0 0 0]),[1,1,1,numel(Kr)]);
        end
        j             = isnan(tmp);
        tmp(j)        = 0;
        msk           = msk + tmp;
    end
    
    uname     = sprintf('U_Sub_%s',num2str(s));
    U         = createvol(uname,vdir,u,'mat',MG,'descrip','Displacemet Field');
    fntu{s}   = U.dat.fname;
end

if all(msk(:)==0)
    msk = 1;
else
    msk(msk>=5)       = 5;
    msk               = msk/5;
    tmp               = sum(msk(:,:,:,1:K1),4);
    tmp(tmp>1)        = 1;
    msk(:,:,:,1:K1-1) = repmat(Smooth(tmp,10,mm),[1,1,1,K1-1]);
    tmp               = 1 - sum(msk(:,:,:,1:K1),4);
    tmp(tmp<0)        = 0;
    msk(:,:,:,K1+1)   = Smooth(tmp,10,mm);
    msk(:,:,:,K1)     = ones(size(msk(:,:,:,K1)));
end


