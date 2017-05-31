
function [y,E]=warp(x,M,phi,prs,varargin)
if ~isempty(varargin)
    d = varargin{1};
    if numel(varargin)>1
        st = varargin{2};
    else
        st = [1 1 1];
    end
    phi = Identity(d,st);
end
E = aff(phi,M);
y = zeros(size(E(:,:,:,1)),'single');
for cl=1:size(x,4)
    x(:,:,:,cl) = spm_diffeo('bsplinc',x(:,:,:,cl),prs);
    y(:,:,:,cl) = spm_diffeo('bsplins',x(:,:,:,cl),E,prs);
end

end