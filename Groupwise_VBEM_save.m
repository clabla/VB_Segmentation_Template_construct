
function res = Groupwise_VBEM_save(In,biaspar,affpar,po,pr,mom,w) 

S   = In.S;
K   = In.K;
C   = size(In.pr.m,1);
prm = In.prm;
prs = [1 1 1 0 0 0];
ist = [1 1 1];
Ker = spm_shoot_greens('kernel',In.DG,prm);

delete(gcp);parpool(3);

parfor s=1:S
    
    iy1 = []; D = []; Mr = []; nr = []; X = []; Xd = []; Bfi = []; di = []; n = []; tpm = [];
    
    B           = nifti(In.fntb{s});
    U           = nifti(In.fntu{s});
    u           = single(U.dat(:,:,:,:));
    TPMs        = nifti(In.TPMs);
    b           = single(TPMs.dat(:,:,:,:));
    [ph,~,~,th] = spm_shoot3d(u,prm,8,Ker);
    
    for c=1:C
        if ~isempty(In.fn{s,c})
            Nii   = nifti(In.fn{s,c});
            [di,n] = fileparts(In.fn{s,c});
            
            if c==1
                nr  = n;
                Mr  = Nii.mat;
                D   = size(Nii.dat(:,:,:));
                y   = Aff(ph,inv(Mtransf(affpar{s},In.MG,eye(4))));
                iy  = Aff(spm_diffeo('comp',th,Aff(Identity(In.DG,ist),Mtransf(affpar{s},In.MG,In.MG))),In.MG);
                y1  = Warp(y,In.MG\Mr,[],[4 4 4 0 0 0],D,ist);
                iy1 = Warp(iy,In.MG\In.MG1,[],[4 4 4 0 0 0],In.DG1,ist);
                u1  = Warp(u,In.MG\In.MG1,[],prs,In.DG1,ist);
                u1  = bsxfun(@times,u1,reshape(sqrt(sum(In.MG(1:3,1:3).^2))./sqrt(sum(In.MG1(1:3,1:3).^2)),[1 1 1 3]));
                JD  = spm_diffeo('def2det',iy1);
                tpm = Warp_template(b,inv(In.MG),y1,prs);
                tpm = Rw_priors(tpm,w{s});                
                Createvol(strcat('u',n),In.resdir,u1,'mat',In.MG1,'descrip','Velocity Field');
            end
            
            if numel(size(Nii.dat))==3
                X  = spm_diffeo('bsplins',single(Nii.dat(:,:,:)),Aff(iy1,inv(Nii.mat)),prs);
                Xd = cat(2,Xd,reshape(Warp(single(Nii.dat(:,:,:)),Nii.mat\Mr,[],prs,D,ist),[prod(D),1]));
            elseif numel(size(Nii.dat))==4
                X  = spm_diffeo('bsplins',single(Nii.dat(:,:,:,c)),Aff(iy1,inv(Nii.mat)),prs);
                Xd = cat(2,Xd,reshape(Warp(single(Nii.dat(:,:,:,c)),Nii.mat\Mr,[],prs,D,ist),[prod(D),1]));
            end
            
            Bf = spm_diffeo('bsplins',single(B.dat(:,:,:,c)),Aff(iy1,inv(B.mat)),prs);
            
            if strcmp(n,nr)
                Createvol(strcat('w',n,'_',num2str(c)),In.resdir,Bf.*X,'mat',In.MG1);
            else
                Createvol(strcat('w',n),In.resdir,Bf.*X,'mat',In.MG1);
            end
            
            Bfi = cat(2,Bfi,reshape(Warp(single(B.dat(:,:,:,c)),B.mat\Mr,[],prs,D,ist),[prod(D),1]));
        else
            Xd = cat(2,Xd,NaN(size(Xd,1),1));
            Bfi = cat(2,Bfi,ones(size(Bfi,1),1));
        end
        
         if exist(fullfile(di,strcat(n,'_labels','.nii')),'file')
            R                = nifti(fullfile(di,strcat(n,'_labels','.nii')));
            wR               = Warp_template(single(R.dat(:,:,:,:)),inv(Mr),iy1,prs*0);
            Createvol(strcat('w',nr,'_labels'),In.resdir,wR,'mat',In.MG1);
        end
        if exist(fullfile(di,strcat(n,'_brain_labels','.nii')),'file')
            R                = nifti(fullfile(di,strcat(n,'_brain_labels','.nii')));
            wR               = Warp_template(single(R.dat(:,:,:,:)),inv(Mr),iy1,prs);
            Createvol(strcat('w',nr,'_brain_labels'),In.resdir,wR,'mat',In.MG1);
        end
         if exist(fullfile(di,strcat(n,'_brain_labels_test','.nii')),'file')
            R                = nifti(fullfile(di,strcat(n,'_brain_labels_test','.nii')));
            wR               = Warp_template(single(R.dat(:,:,:,:)),inv(Mr),iy1,prs);
            Createvol(strcat('w',nr,'_brain_labels_test'),In.resdir,wR,'mat',In.MG1);
        end
        if exist(fullfile(di,strcat(n,'_spine_labels','.nii')),'file')
            R                = nifti(fullfile(di,strcat(n,'_spine_labels','.nii')));
            wR               = Warp_template(single(R.dat(:,:,:,:)),inv(Mr),iy1,prs);
            Createvol(strcat('w',nr,'_spine_labels'),In.resdir,wR,'mat',In.MG1);
        end
        if exist(fullfile(di,strcat(n,'_spine_labels_test','.nii')),'file')
            R                = nifti(fullfile(di,strcat(n,'_spine_labels_test','.nii')));
            wR               = Warp_template(single(R.dat(:,:,:,:)),inv(Mr),iy1,prs);
            Createvol(strcat('w',nr,'_spine_labels_test'),In.resdir,wR,'mat',In.MG1);
        end
        
    end
    
    rs = rRsp(Xd,Bfi,reshape(tpm,[prod(D),K]),po{s});
    rs(:,4:end) = [];
    rs(sum(rs(:,1:3),2)<=0.75,1:3) = 0;
    rs = reshape(rs,[D size(rs,2)]);
    wR = Warp_template(rs,inv(Mr),iy1,prs);
    wR = bsxfun(@times,wR,JD);
    
        Createvol(strcat('seg',nr),In.resdir,rs,'mat',Mr);
        Createvol(strcat('wpr',nr),In.resdir,tpm,'mat',Mr);
        Createvol(nr,In.resdir,reshape(Xd.*Bfi,[D,size(Xd,2)]),'mat',Mr);
        Createvol(strcat('segw',nr),In.resdir,wR,'mat',In.MG1);
    
end
res = struct('posteriors',reshape(po,S,1),...
             'priors'    ,pr,...
             'bias'      ,reshape(biaspar,S,1),...
             'affine'    ,reshape(affpar,S,1),...
             'mom'       ,reshape(mom,S,1),...
             'w'         ,reshape(w,S,1));
save(fullfile(In.resdir,'res'),'res','-v7.3')
n = fullfile(In.resdir,'input_files');
save(n,'In','-v7.3')
end

