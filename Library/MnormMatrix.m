function [x] = MnormMatrix(M, A, N);
    
S = M*A*N;
x = sqrt(svds(S' * S, 1));
% slows down the calculation:
% check = abs(x-norm(full(M*A*N)));
% if check > 10^-15
%     warning('difference between norm and svds is    : %e \n\n',check)
% end
    
if imag(x) > 10^-15
    warning('Imaginray part to big')
    imag(x)
end
x = real(x);
end