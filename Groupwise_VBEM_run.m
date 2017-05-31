% C  -> number of channels;
% K  -> number of clusters;
% Brain spine merged

function [In,affpar,biaspar,po,pr,mom,w] = GW_diff_VBEM_multlab_par_run(In)

K       = In.K;
Km      = find(In.tiss_ind == In.K1,1,'last');
C       = In.C;
S       = In.S;
a       = In.a;
am      = In.a;
po      = In.po;
mom     = cell(S,1);
step    = cell(S,1);
l_b     = cell(S,1);
l_d     = cell(S,1);
pr      = cell(S,1);
affpar  = In.affpar;
biaspar = In.biaspar;
w       = In.w;
ist     = [1 1 1];
reg_st  = 30;
L       = 0;
mm      = sqrt(sum(In.MG(1:3,1:3).^2));
prm     = In.prm;
Ker     = spm_shoot_greens('kernel',In.DG,prm);
or      = size(In.biaspar{1}(:,:,:,1));
BN      = prod(or);
N       = prod(In.DG);
for s=1:S
    pr{s}   = In.pr;
    l_d{s}  = 1e-1;
    l_b{s}  = 1e+1;
    step{s} = 1;
end
l_dmax = 1e+10;
l_dmin = 1e-10;
l_bmin = 1e+1;
l_bmax = 1e+10;
for s=1:S
    step{s}= 1;                                                           
end
delete(gcp); parpool(6);

%start iterating===========================================================

for iter=1:30
    
    fprintf('Iteration #%s \n',num2str(iter));
    ite1 = max(iter,reg_st+1);
    prm1 = [prm(1:4) prm(5:end)/(1*(1-exp(-((ite1-reg_st)^2)/100)))];
    oL   = L;
    L    = 0;
    msk  = zeros([In.DG,In.K1]);
    
    
    %create temporary file for TPMs----------------------------------------
    
    btmp = zeros([In.DG,K],'single');
    dtmp = zeros([In.DG,K],'single');
    
    for s=1:S %loop over subjects
        
        fprintf('Subject #%s \n',num2str(s));
        
        %load image, bias, In.TPMs and displacement field------------------
        
        I      = nifti(In.fnt{s});
        B      = nifti(In.fntb{s});
        U      = nifti(In.fntu{s});
        TPMs   = nifti(In.TPMs);
        Xm     = reshape(single(I.dat(:,:,:,:)),[N,C]);
        Bf     = reshape(single(B.dat(:,:,:,:)),[N,C]);
        u      = single(U.dat(:,:,:,:)); 
        b      = single(TPMs.dat(:,:,:,:)); 
        if ~isempty(In.fnr{s})
            R  = nifti(In.fnr{s}.name);
            rl = reshape(single(R.dat(:,:,:,:)),[N,size(R.dat,4)]);
            Kr = In.fnr{s}.K;
        else
            rl = [];
            Kr = [];
        end
        Sigma = zeros(C,C,K);
        tmp1  = zeros([In.DG,In.K1]);
        
        %compute deformations and Jacobians from v0------------------------
        
        [ph,Jph]  = spm_shoot3d(u,prm,In.t_arg,Ker);
        Xsi       = Aff(ph,inv(Mtransf(affpar{s},In.MG,In.MG)));
        
        %transform TPMs----------------------------------------------------
        
        bs = Warp_template(b,eye(4),Xsi,In.prs) ;  Xsi = [];  bs = Rw_priors(bs,w{s});
        
        %compute terms for regularization----------------------------------
        
        bp    = reshape(biaspar{s},BN,C);
        Lregd = - 0.5*affpar{s}'*In.Cd^-1*affpar{s} - sum(log(diag(chol(In.Cd)))) - 0.5*C*log(2*pi);
        Lregv = - 0.5*sum(sum(sum(sum(u.*spm_diffeo('vel2mom',u,prm1)))));
        Lregb = - 0.5*bp(:)'*In.Pb*bp(:) - sum(log(diag(chol(In.Pb^-1)))) - 0.5*C*log(2*pi) + sum(log(prod(Bf,2)));
        
        [po{s},pr{s},mom{s},Lk,Elpz,Eqz,r,w{s}] = VB_EM(Xm,Bf,bs,po{s},pr{s},w{s},K,N,In.tiss_ind,rl,Kr);
        
        if 1
            for i=1:2
                
                %compute suff statistics and update MoG hyperparameters----
                [po{s},pr{s},mom{s},Lk,Elpz,Eqz,r,w{s}] = VB_EM(Xm,Bf,bs,po{s},pr{s},w{s},K,N,In.tiss_ind,rl,Kr);
                
            end
        end
        
        lb = nansum(Lk(1,:) + Lk(2,:) - Lk(3,:) - Lk(4,:)) + Elpz - Eqz + Lregb + Lregd + Lregv;
        
        if iter>0
            %estimate bias=================================================
            
            for k=1:K
                Sigma(:,:,k) = inv(po{s}.W(:,:,k))/(po{s}.nu(k)-C-1);
            end
            
            for ib=1:1
                
                olb = lb;
                opo = po{s};
                obp = bp;                        %old bias field parameters
                
                for c=1:C                               %loop over channels                    
                    if ~isempty(In.fn{s,c})
                        l  = zeros(1,C); l(1,c)= 1; l = logical(l);
                        
                        %compute objective function derivatives------------
                        
                        [gb,Hb] = Bias_derivs(In.DG,Xm(:,c),r,po{s}.m(l,:),Sigma(l,l,:),Bf(:,c),or);
                        Cbch    = In.Pb(BN*(c - 1)+1:BN*c,BN*(c - 1)+1:BN*c)*1;
                        gb      = gb + Cbch*obp(:,c);
                        Hb      = Hb + Cbch;
                        Hbl     = l_b{s}*eye(BN);
                        
                        %bias field parameters update----------------------
                        
                        bp(:,c) = obp(:,c) - (Hb + Hbl)\gb;
                    end
                end
                
                %recompute bias--------------------------------------------
                
                Bf = Bias_Field(In.DG,C,or,bp,a{s});
                if any(~isfinite(Bf(:))); warning('Bias not finite'); end
  
                
                %update suff statistics and MoG hyperparameters------------
                
                [po{s},pr{s},mom{s},Lk,Elpz,Eqz,r,w{s}]  = VB_EM(Xm,Bf,bs,po{s},pr{s},w{s},K,N,In.tiss_ind,rl,Kr);
                               
                %compute objective function--------------------------------
                
                Lregb = - 0.5*bp(:)'*In.Pb*bp(:) - sum(log(diag(chol(In.Pb^-1)))) - 0.5*C*log(2*pi) + sum(log(prod(Bf,2)));
                lb    = nansum(Lk(1,:) + Lk(2,:) - Lk(3,:) - Lk(4,:)) + Elpz - Eqz + Lregb + Lregd + Lregv;
                
                if lb>=olb
                    l_b{s} = max(l_b{s}/10,l_bmin);
                    break
                elseif lb<olb
                    l_b{s}                       = min(l_b{s}*10,l_bmax);
                    lb                           = olb;
                    bp                           = obp;
                    po{s}                        = opo;
                    Bf                           = Bias_Field(In.DG,C,or,bp,a{s});
                    [po{s},pr{s},mom{s},~,~,~,r] = VB_EM(Xm,Bf,bs,po{s},pr{s},w{s},K,N,In.tiss_ind,rl,Kr);
                end
            end
            
            am{s}          = nanmean(Bf,1);
            biaspar{s}     = reshape(bp,[or,C]);
            B.dat(:,:,:,:) = reshape(Bf,[In.DG,C]);
        end
        
        if iter>2
            %estimate affine deformations==================================
            
            a_s = affpar{s};
            
            for id=1:2
                
                olb  = lb;
                oa_s = a_s;
                
                %compute log-likelihood derivatives------------------------
                
                [M,dM]  = Mtransf(a_s,In.MG,In.MG);
                [gd,Hd] = AffDerivs(r,b,inv(M),dM,ph,In.Cd,oa_s,w{s});
                
                %deformation parameters update-----------------------------
                
                a_s = oa_s - (Hd + l_d{s}*eye(12))\gd;
                
                %transform In.TPMs-----------------------------------------
                Xsi = Aff(ph,inv(Mtransf(a_s,In.MG,In.MG)));
                bs  = Warp_template(b,eye(4),Xsi,In.prs);  Xsi = [];  bs = Rw_priors(bs,w{s});
                
                %update suff statistics and MoG hyperparameters------------
                
                [po{s},pr{s},mom{s},Lk,Elpz,Eqz,r,w{s}]  = VB_EM(Xm,Bf,bs,po{s},pr{s},w{s},K,N,In.tiss_ind,rl,Kr);
                
                
                %compute objective function--------------------------------
                
                Lregd = - 0.5*a_s'*In.Cd^-1*a_s - sum(log(diag(chol(In.Cd)))) - 0.5*C*log(2*pi);
                lb    = nansum(Lk(1,:) + Lk(2,:) - Lk(3,:) - Lk(4,:)) + Elpz - Eqz + Lregb + Lregd + Lregv;
                
                if lb>=olb
                    l_d{s} = max(l_dmin,l_d{s}/10);
                    break
                else
                    l_d{s}                       = min(l_d{s}*10,l_dmax);
                    a_s                          = oa_s;
                    lb                           = olb;
                    Xsi                          = Aff(ph,inv(Mtransf(a_s,In.MG,In.MG)));
                    bs                           = Warp_template(b,eye(4),Xsi,In.prs);  bs = Rw_priors(bs,w{s});
                    [po{s},pr{s},mom{s},~,~,~,r] = VB_EM(Xm,Bf,bs,po{s},pr{s},w{s},K,N,In.tiss_ind,rl,Kr);
                end
            end
            
            affpar{s} = a_s;
        end
        
        if iter>reg_st
            %estimate nonlinear deformations===============================
            
            M = Mtransf(affpar{s},In.MG,In.MG);
            [gd,Hd]  = VDerivs(In.DG,r,b,w{s},inv(M),ph,Jph); 
            gd       = gd + spm_diffeo('vel2mom',u,prm1);
            dv       = spm_diffeo('fmg',Hd,gd,[prm1,3,2]);
            
            for id=1:2
                
                olb    = lb;
                opo    = po{s};
                oLregv = Lregv;
                u      = u - step{s}*dv;  
                ph     = spm_shoot3d(u,prm,In.t_arg,Ker);
                Xsi    = Aff(ph,inv(Mtransf(affpar{s},In.MG,In.MG)));
                
                %transform TPMs--------------------------------------------
                
                bs = Warp_template(b,eye(4),Xsi,In.prs);  bs = Rw_priors(bs,w{s});
                
                %update suff statistics and MoG hyperparameters------------
                
                [po{s},pr{s},mom{s},Lk,Elpz,Eqz,r,w{s}]  = VB_EM(Xm,Bf,bs,po{s},pr{s},w{s},K,N,In.tiss_ind,rl,Kr);
                
                %compute objective function--------------------------------
                
                Lregv = - 0.5*sum(sum(sum(sum(u.*spm_diffeo('vel2mom',u,prm1)))));
                lb = nansum(Lk(1,:) + Lk(2,:) - Lk(3,:) - Lk(4,:)) + Elpz - Eqz + Lregb + Lregd + Lregv;
                
                if 1
                    if lb>=olb
                        step{s} = min(step{s}*1.1,1);
                        break
                    else
                        u  = single(U.dat(:,:,:,:));
                        Lregv = oLregv;
                        step{s}                                 = max(step{s}/1.1,0.2);
                        po{s}                                   = opo;
                        ph                                      = spm_shoot3d(u,prm,In.t_arg,Ker);
                        Xsi                                     = Aff(ph,inv(Mtransf(affpar{s},In.MG,In.MG)));
                        bs                                      = Warp_template(b,eye(4),Xsi,In.prs);  
                        bs                                      = Rw_priors(bs,w{s});
                        [po{s},pr{s},mom{s},Lk,Elpz,Eqz,r,w{s}] = VB_EM(Xm,Bf,bs,po{s},pr{s},w{s},K,N,In.tiss_ind,rl,Kr);           
                        lb                                      = nansum(Lk(1,:) + Lk(2,:) - Lk(3,:) - Lk(4,:)) + Elpz - Eqz + Lregb + Lregd + Lregv;
                    end
                end
            end
            U.dat(:,:,:,:) = u;
        end
        
        if In.tpm_upd 
            if iter>0
            %update statistics for TPMs====================================
            
            [~,~,~,ph,Jph]   = spm_shoot3d(u,prm,In.t_arg,Ker);
            ph               = spm_diffeo('comp',ph,Aff(Identity(In.DG,ist),Mtransf(affpar{s},In.MG,In.MG)));
            tmp              = Warp_template(single(reshape(r,[In.DG,K])),eye(4),ph,In.prs);
            tmp(isnan(tmp))  = 0;
                       
            if  ~isempty(In.fnr{s})
                if numel(In.fnr{s}.K)==size(rl,2)
                    tmp1(:,:,:,In.fnr{s}.K) = Warp_template(reshape(rl,[In.DG,size(rl,2)]),eye(4),ph,In.prs);
                elseif size(rl,2)==1 && numel(In.fnr{s}.K)>1
                    tmp1(:,:,:,In.fnr{s}.K) = repmat(Warp_template(reshape(rl,[In.DG,size(rl,2)]),eye(4),ph,In.prs),[1,1,1,numel(In.fnr{s}.K)]);
                end
                j       = isnan(tmp1);
                tmp1(j) = 0;
                msk     = msk + tmp1;
            end
            tmp1 = [];
            
            if numel(size(In.msk)) == numel(size(tmp(:,:,:,1:In.K1)))
                if size(In.msk,4)-1  == size(tmp(:,:,:,1:In.K1),4)
                    tmp(:,:,:,1:Km)   = tmp(:,:,:,1:Km).*In.msk(:,:,:,In.tiss_ind(1:Km));
                    tmp(:,:,:,Km+1:K) = bsxfun(@times,tmp(:,:,:,Km+1:K),In.msk(:,:,:,end));
                end
            end
            
            
            Jd      = spm_diffeo('det',Jph)*det(Mtransf(affpar{s},In.MG,In.MG));
            tmp1    = bsxfun(@times,tmp,Jd) + In.alpha0;
            j       = isnan(tmp1);
            tmp1(j) = 0;
            btmp    = btmp + tmp1;
            tmp1    = bsxfun(@times,Rw_priors(b,w{s})./b,sum(tmp,4).*Jd) + In.alpha0;
            j       = isnan(tmp1);
            tmp1(j) = 0;
            dtmp    = dtmp + tmp1;
            tmp     = [];
            tmp1    = [];
            
            
            %update global objective function==============================
            
            L = L + lb;
            end
        end
        
    end
    
    amm = prod(cell2mat(am),1).^(1/size(am,1));
    for s = 1:S
        a{s} = a{s}./amm;
    end
    
    %normalize In.TPMs-----------------------------------------------------
    
    if In.tpm_upd
        if iter>0
            TPMs                  = nifti(In.TPMs);
            b                     = bsxfun(@rdivide,btmp,dtmp + eps);
            b                     = bsxfun(@rdivide,b,sum(b,4) + eps);
            TPMs.dat(:,:,:,:)     = b;
            if all(msk(:)==0)
                msk = 1;
            else
                msk(msk>=5)           = 5;
                msk                   = msk/5;
                for k=1:In.K1
                    msk(:,:,:,k) = smooth(msk(:,:,:,k),10,mm);
                end
                tmp                   = 1 - sum(msk(:,:,:,1:In.K1),4);
                tmp(tmp<0)            = 0;
                msk(:,:,:,In.K1+1)    = smooth(tmp,10,mm);
                msk(:,:,:,In.K1)      = ones(size(msk(:,:,:,In.K1)));
            end
            In.msk = msk;
        end
    end
    
    %compute intensity priors----------------------------------------------
    for s=1:S
        for c=1:C
            fn(s,c) =~ isempty(In.fn{s,c});
        end
    end
    prn = priors(po,pr{1},find(sum(fn,2)==C));
    for s=1:S
        pr{s} = prn;
    end
    
    %check for convergence-------------------------------------------------
    
    if iter>=2
        if abs((L - oL)/oL)<=1e-6
            break
        end
    end
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