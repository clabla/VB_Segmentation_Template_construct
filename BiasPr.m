%generate precision matrix for bias field parameters-----------------------

function C=BiasPr_v2(dim,or,mm)
Bx = spm_dctmtx(dim(1),or(1));
By = spm_dctmtx(dim(2),or(2));
Bz = spm_dctmtx(dim(3),or(3));

D2x          = sparse(toeplitz([2 -1 zeros(1,dim(1)-3) -1])); 
D2x(1,end)   = 0;  
D2x(end,1)   = 0; 
D2x          = D2x/mm(1)^2;

D2y          = sparse(toeplitz([2 -1 zeros(1,dim(2)-3) -1])); 
D2y(1,end)   = 0; 
D2y(end,1)   = 0; 
D2y          = D2y/mm(2)^2;

D2z          = sparse(toeplitz([2 -1 zeros(1,dim(3)-3) -1])); 
D2z(1,end)   = 0; 
D2z(end,1)   = 0; 
D2z          = D2z/mm(3)^2;

C   = kron((D2z*Bz)'*(D2z*Bz),kron(By'*By,Bx'*Bx))+...
      kron(Bz'*Bz,kron((D2y*By)'*(D2y*By),Bx'*Bx))+...
      kron(Bz'*Bz,kron(By'*By,(D2x*Bx)'*(D2x*Bx)))+...
      2*kron((D2z*Bz)'*Bz,kron((D2y*By)'*By,Bx'*Bx))+...
      2*kron((D2z*Bz)'*Bz,kron(By'*By,(D2x*Bx)'*Bx))+...
      2*kron(Bz'*Bz,kron((D2y*By)'*By,(D2x*Bx)'*Bx));
  C = C/1e+2;
  C(C==0) = eps*eps;
end