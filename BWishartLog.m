
%Logarithm of normalization constant for Wishart Distribution

function logB=BWishartLog(W,nu)

D         = size(W,1);
prodgamln = sum(gammaln(0.5*(nu + 1 - [1:D])));
logB      = - 0.5*nu*logdet(W) - 0.5*nu*D*log(2) - 0.25*D*(D-1)*log(pi) - prodgamln;