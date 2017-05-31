
function Bf = Bias_Field(D,C,or,bp,a)

Bf = zeros([D C]);
Bx = spm_dctmtx(D(1),or(1));
By = spm_dctmtx(D(2),or(2));
Bz = spm_dctmtx(D(3),or(3));

for c=1:C
    for z = 1:D(3)
        tmp = Bx*sum(bsxfun(@times,permute(Bz(z,:),[3 1 2]),reshape(bp(:,c),or)),3)*By';
        Bf(:,:,z,c) = exp(tmp)*a(c);
    end
end
Bf = reshape(Bf,[prod(D),C]);