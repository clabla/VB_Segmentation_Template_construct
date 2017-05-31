function [fnt, fnr] = Init_im(Kbrain,DG,MG,C,S,fn,use_lab)

fnt = cell(S,1);
fnr = cell(S,1);

parfor s=1:S %loop over subjects
    
    K1        = Kbrain;
    Xm        = NaN([DG,C]);
    di        = [];
    na        = [];
    
    for c=1:C %loop over channels
        if ~isempty(fn{s,c})
            Nii     = nifti(fn{s,c});
            [di,na] = fileparts(fn{s,c});
            M = Nii.mat;
            if size(size(Nii.dat),2)==3
                X           = single(Nii.dat(:,:,:));
                Xm(:,:,:,c) = Warp(X,M\MG,[],[0 0 0 0 0 0],DG,[1 1 1]);
            elseif size(size(Nii.dat),2)==4
                X           = single(Nii.dat(:,:,:,c));
                Xm(:,:,:,c) = Warp(X,M\MG,[],[0 0 0 0 0 0],DG,[1 1 1]);
            end
        end
        
        if use_lab
            if exist(fullfile(di,strcat(na,'_brain_labels','.nii')),'file') && exist(fullfile(di,strcat(na,'_spine_labels','.nii')),'file')
                R                = nifti(fullfile(di,strcat(na,'_brain_labels','.nii')));
                rs               = Warp(single(R.dat(:,:,:,:)),R.mat\MG,Identity(DG,[1 1 1]), [0,0,0,0,0,0]);
                R                = nifti(fullfile(di,strcat(na,'_spine_labels','.nii')));
                rs(:,:,:,1:K1-1) = bsxfun(@plus,rs(:,:,:,1:K1-1),Warp(single(R.dat(:,:,:,:)),R.mat\MG,Identity(DG,[1 1 1]),[0,0,0,0,0,0]));
                rdir             = fullfile(tempdir,'\Work\R');
                rname            = sprintf('R_Sub_%s',num2str(s));
                R                = Createvol(rname,rdir,rs,'mat',MG);
                fnr{s}.name      = R.dat.fname;
                fnr{s}.labeltype = 'brain&spine';
                fnr{s}.K = 1:K1;
            elseif exist(fullfile(di,strcat(na,'_brain_labels','.nii')),'file')
                R                = nifti(fullfile(di,strcat(na,'_brain_labels','.nii')));
                rs               = Warp(single(R.dat(:,:,:,:)),R.mat\MG,Identity(DG,[1 1 1]),[0,0,0,0,0,0]);
                rdir             = fullfile(tempdir,'\Work\R');
                rname            = sprintf('R_Sub_%s',num2str(s));
                R                = Createvol(rname,rdir,bsxfun(@rdivide,rs,sum(rs,4)+eps),'mat',MG);
                fnr{s}.name      = R.dat.fname;
                fnr{s}.labeltype = 'brain';
                fnr{s}.K = 1:K1;
                
            elseif exist(fullfile(di,strcat(na,'_spine_labels','.nii')),'file')
                R                = nifti(fullfile(di,strcat(na,'_spine_labels','.nii')));
                rs               = Warp(single(R.dat(:,:,:,:)),R.mat\MG,Identity(DG,[1 1 1]),[0,0,0,0,0,0]);
                rdir             = fullfile(tempdir,'\Work\R');
                rname            = sprintf('R_Sub_%s',num2str(s));
                R                = Createvol(rname,rdir,repmat(rs,[1,1,1,K1-1]),'mat',MG);
                fnr{s}.name      = R.dat.fname;
                fnr{s}.labeltype = 'spine';
                fnr{s}.K = 1:K1-1;
            end
        end
    end
    
    imdir = fullfile(tempdir,'\Work\Im');
    if ~isdir(imdir)
        mkdir(imdir)
    end
    
    imname = sprintf('Im_Sub_%s',num2str(s));
    I      = Createvol(imname,imdir,Xm,'mat',MG);
    fnt{s} = I.dat.fname;
end