function In = Groupwise_VBEM_initialize(In) 

%define useful variables and parameters====================================

K         = In.K;
Kbrain    = 3;
C         = In.C;
S         = In.S;
mm        = [3 3 3];                                                          
ist       = [1 1 1];                                                    
In.prs    = [1,1,1,0,0,0];                                                    
In.t_arg  = 8;
a         = cell(S,1);
w         = cell(S,1);
In.l_d    = cell(S,1);
In.l_b    = cell(S,1);                                                         
In.Cd(1:3,1:3) = 1e-0*eye(3); In.Cd(4:6,4:6) = 1e-3*eye(3); In.Cd(7:9,7:9) = 0.5e-3*eye(3); In.Cd(10:12,10:12) = 0.5e-5*eye(3);                                                   
                                                  
resdir = fullfile(In.di,sprintf('%ssub%sclass%s_%s_%smm',num2str(S),num2str(K),num2str(mm(1)),num2str(mm(2)),num2str(mm(3))));
if ~isdir(resdir)
    mkdir(resdir);
end
In.resdir = resdir;

if ~isfield(In,'use_lab')
    In.use_lab = 0;
end

%reslice channels==========================================================

MF  = [];
DF  = [];

for s=1:S
    ch  = [];
    for c=1:C
        ch = logical([ch,~isempty(In.fn{s,c})]);
    end
    Nii = nifti(In.fn(s,ch));
    for c=1:sum(ch)
        MF = cat(3,MF,Nii(c).mat);
        DF = cat(1,DF,size(Nii(c).dat));
    end
end

if ~isfield(In,'TPMs0')
    [MG1,DG1] = Compute_avg_mat(MF,DF); clear MF DF
else
    TPMs0     = nifti(In.TPMs0);
    MG1       = TPMs0.mat;
    DG1       = TPMs0.dat.dim(1:3);
    [MG1,DG1] = Compute_avg_mat(cat(3,MG1,MF),cat(1,DG1,DF)); clear MF DF
end

vx        = sqrt(sum(MG1(1:3,1:3).^2));
st        = mm./vx;
F         = diag(st);
F(1:4,4)  = 1;
MG        = MG1*F;
if isfield(In,'udir')
    DG = floor(DG1./st) - 1;
else
    DG = floor(DG1./st);
end
N         = prod(DG);

[fnt, fnr] = Init_im(Kbrain,DG,MG,C,S,In.fn,In.use_lab);

%initialize intensity distribution hyperparameters=========================

if ~isfield(In,'a')
    for s=1:S
        I    = nifti(fnt(s));
        for c=1:C
            h = I.dat(:,:,:,c);
            h(h<=1) = NaN;
            m = nanmean(reshape(h,[N,1]));
            if isnan(m); m = 512; end;
            a{s}(c) = 512/m;  
        end
        
    end
else
    a = In.a;
end
if isfield(In,'pr0')
    pr = In.pr0;
else
    init{1}  = 'rand';
    init{2}  = 150;
    s = randsample(1:S,1);
    for c=1:C;
        while isempty(In.fn{s,c})
            s = randsample(1:S,1);
        end
        I   = nifti(fnt(s));
        Xm(:,c)  = bsxfun(@times,reshape(single(I.dat(:,:,:,c)),[N,1]),a{s}(c));        
    end
    [mu,Sig] = Kmeans(Xm,K,init);
    pr.beta  = 0.01*ones(1,K);
    pr.m     = mu;
    pr.nu    = (C)*ones(1,K)-0.99;
    for k=1:K
        pr.W(:,:,k) = inv(Sig(:,:,k))/pr.nu(k)/1e+1;
    end
end
 po = cell(S,1);
if isfield(In,'po')
    if ~isempty(In.po{1}) && size(In.po{1}.beta,2) == K 
        po = In.po;
    else
        if ~isfield(In,'TPMs0')          
            for s=1:S
                po{s} = pr;
            end
        end
    end
elseif ~isfield(In,'po') && ~isfield(In,'TPMs0')
    for s=1:S
        po{s} = pr;
    end
end
if ~isfield(In,'w')
    w = cell(S,1);
    for s=1:S
        w{s} = ones(1,K);
    end
end

%initialize TPMS===========================================================

alpha0 = 0.01;
if isfield(In,'TPMs0')
    TPMs0  = nifti(In.TPMs0);
    b      = Warp(single(TPMs0.dat(:,:,:,:)),TPMs0.mat\MG,[],[1 1 1 0 0 0],DG,ist);
    b(b<0) = 0;
    if size(b,4)<K
        b(:,:,:,end+1:K) = 1/K;
        b = bsxfun(@rdivide,b,sum(b,4));
    elseif size(b,4)>K
        b(:,:,:,K+1:end) = [];
        b = bsxfun(@rdivide,b,sum(b,4));
    else
        b = bsxfun(@rdivide,b+eps,sum(b+eps,4));
    end
else
    b = ones(prod(DG),K,'single');
    b = reshape(bsxfun(@rdivide,b,sum(b,2)),[DG,K]);   
end
tpmdir  = fullfile(tempdir,'Work\tpms');
if ~isdir(tpmdir); mkdir(tpmdir); end
Createvol('tpms',resdir,b(:,:,:,1:K),'mat',MG,'descrip','tpms');

%initialize bias field=====================================================

fntb    = cell(S,1);
In.biaspar = cell(S,1);
if isfield(In,'bias')
    In.biaspar  = In.bias;
    or = size(In.bias{1}(:,:,:,1));
else
    or = [5 5 5];   
    for s=1:S
    In.biaspar{s}  = zeros([or,C]);   
    end
end

BN    = prod(or);
bfdir = fullfile(tempdir,'\Work\Bf');
if ~isdir(bfdir); mkdir(bfdir); end

for c=1:C
        In.Pb((c-1)*BN+1:c*BN,(c-1)*BN+1:c*BN) = BiasPr_v2(DG,or,mm);            
end

parfor s=1:S  
    Bf      = Bias_Field(DG,C,or,reshape(In.biaspar{s},BN,C),a{s});
    bfname  = sprintf('Bf_Sub_%s',num2str(s));
    B       = Createvol(bfname,bfdir,reshape(Bf,[DG,C]),'mat0',MG,'descrip','Bias Field');
    fntb{s} = B.dat.fname;
end

clear B Bf Bx By Bz


%initialize displacements==================================================
 
mm      = sqrt(sum(MG(1:3,1:3).^2));
In.prm  = [mm 0.00005 0.8*[0.005 1 0.2 1]];

vdir = fullfile(tempdir,'\Work\V');

if ~isdir(vdir); mkdir(vdir); end

In.affpar = cell(S,1);

if isfield(In,'affine')
    In.affpar = In.affine;
else
    for s=1:S
        In.affpar{s}=zeros(12,1);
    end
end

if isfield(In,'udir')
    udir = In.udir;
else
    udir = [];
end

[fntu, msk] = Init_u(DG,MG,C,Kbrain,S,In.fn,fnr,udir,vdir,mm,In.affpar,In.prm);

In.ngaus    = ones(1,Kbrain);
In.ngaus(1) = 2;

if (Kbrain)>1
    In.tiss_ind = ones(1,In.ngaus(1));
    for k=2: Kbrain
        In.tiss_ind = cat(2, In.tiss_ind ,repmat(k,[1,In.ngaus(k)]));
    end
    for k=numel(In.tiss_ind)+1:In.K
        In.tiss_ind(k) =  In.tiss_ind(k-1) + 1;
    end
else
    In.tiss_ind = 1:In.K;
end

if ~isfield(In,'tpm_upd')
    In.tpm_upd = 1;
end
In.DG     = DG;
In.DG1    = DG1;
In.MG     = MG;
In.MG1    = MG1;
In.K1     = Kbrain;
In.fnt    = fnt;
In.fntb   = fntb;
In.fntu   = fntu;
In.fnr    = fnr;
In.w      = w;
In.a      = a;
In.po     = po;
In.pr     = pr;
In.alpha0 = alpha0;
In.msk    = msk;
In.TPMs   = fullfile(resdir,'tpms.nii');
end