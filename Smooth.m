function out=smooth(in,fwhm,mm)

out  = zeros(size(in));
fwhm =  fwhm/nthroot(prod(mm),3);
il   = round(3*fwhm);
h    = spm_smoothkern(fwhm,-il:il,1);
spm_conv_vol(single(in),out,h,h,h,-[il il il]);

end
