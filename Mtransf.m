function [M,dMinv,dM]=Mtransf(par,MF,MG)

dM    = zeros(12,12);
dMinv = zeros(12,12);

logT  = par;
dlogT = zeros(4,4,numel(logT));
dlogT(1,4,1)  = 1;
dlogT(2,4,2)  = 1;
dlogT(3,4,3)  = 1;
dlogT(2,3,4)  = 1; dlogT(3,2,4) = -1;
dlogT(3,1,5)  = 1; dlogT(1,3,5) = -1;
dlogT(1,2,6)  = 1; dlogT(2,1,6) = -1;
dlogT(1,1,7)  = 1;
dlogT(2,2,8)  = 1;
dlogT(3,3,9)  = 1;
dlogT(1,2,10) = 1; dlogT(2,1,10)  = 1;
dlogT(1,3,11) = 1; dlogT(3,1,11)  = 1;
dlogT(2,3,12) = 1; dlogT(3,2,12)  = 1;

[T,dTtmp] = spm_dexpm(logT,dlogT);
M = MF\T*MG;

for g=1:size(dTtmp,3)
    dMtmp      = MF\dTtmp(:,:,g)*MG;
    dMinvtmp   = -M\dMtmp/M;
    dMtmp      = dMtmp';
    dMinvtmp   = dMinvtmp';
    dMtmp      = dMtmp(:,1:3);
    dMinvtmp   = dMinvtmp(:,1:3);
    dM(:,g)    = dMtmp(:);
    dMinv(:,g) = dMinvtmp(:);
end
end