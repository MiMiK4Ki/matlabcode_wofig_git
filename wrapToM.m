function lambda = wrapToM(lambda,M)
% wrapToM
%   実数は従来の wrap を適用。
%   複素数は実部・虚部を独立に wrap する。

if ~isreal(lambda)
    lambda = wrapToM(real(lambda), M) + 1j * wrapToM(imag(lambda), M);
    return;
end

lambda = lambda + M;

positiveInput = (lambda > 0);
lambda = mod(lambda, 2*M);
lambda((lambda == 0) & positiveInput) = 2*M;

lambda = lambda - M;

end
