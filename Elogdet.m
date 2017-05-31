
%Expectation of logarithm of determinant of Wishart distributed covariance matrix

function lLambda = Elogdet(W,nu)

C       = size(W,1);
ldW     = 2*sum(log(diag(chol(W))));
Sdig    = sum(psi((1 + nu-(1:C))/2));
lLambda = ldW + C*log(2) + Sdig;