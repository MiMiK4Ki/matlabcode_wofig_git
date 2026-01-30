function lambda = wrapToM(lambda,M)

lambda = lambda +M;

positiveInput = (lambda > 0);
lambda = mod(lambda, 2*M);
lambda((lambda == 0) & positiveInput) = 2*M;

lambda = lambda - M;

end