
function l = logdet(W)

l = 2*sum(log(diag(chol(W))));