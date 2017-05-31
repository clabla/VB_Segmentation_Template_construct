
function bw = Rw_priors(b,w)

K = size(b,4);
N = numel(b(:,:,:,1));
if isempty(w)
    w = ones(1,K);
elseif size(w,1)==K
    w  = reshape(w,1,K);
end
bw = bsxfun(@times,reshape(b,[N,K]),w);
bw = reshape(bw,size(b));
bw = bsxfun(@rdivide,bw,sum(bw,4)+eps);