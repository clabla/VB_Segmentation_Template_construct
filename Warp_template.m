function [y,th]=warp_template(x,M,th,prs)

if ~isempty(M)
    th= aff(th,M);
end
y = zeros(size(th(:,:,:,1)),'single');
for cl=1:size(x,4)
    y(:,:,:,cl) = spm_diffeo('bsplins',x(:,:,:,cl),th,prs);
end
end